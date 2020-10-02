from firedrake import *
import numpy
import subprocess
import os
import warnings
from firedrake.mg.interface import prolong
import numpy as np


def render(mesh,i):
    V = FunctionSpace(mesh, "CG", 1)
    u = Function(V)
    File("test/level%d.pvd"%i).write(u)
    return

def create_labeled_hierarchy(coarse, stepfile,levels, order=1, comm=COMM_WORLD, distribution_parameters=None, callbacks=None, reorder = None):
    try:
        from OCC.Core.STEPControl import STEPControl_Reader
        from OCC.Extend.TopologyUtils import TopologyExplorer
    except ImportError:
        raise ImportError("To use OpenCascadeMeshHierarchy, you must install firedrake with the OpenCascade python bindings (firedrake-update --opencascade).")

    if not os.path.isfile(stepfile):
        raise OSError("%s does not exist" % stepfile)

    project_refinements_to_cad = True

    step_reader = STEPControl_Reader()
    step_reader.ReadFile(stepfile)
    step_reader.TransferRoot()
    shape = step_reader.Shape()
    cad = TopologyExplorer(shape)

    mh = MeshHierarchy(coarse,levels, distribution_parameters=distribution_parameters, callbacks=callbacks, reorder=reorder)
    
    project_to_cad = project_mesh_to_cad_2d
    
    for mesh in mh:
        project_to_cad(mesh, cad)
    mh.nested = False

    if order > 1:
        VFS = VectorFunctionSpace
        Ts = [
            Function(VFS(mesh, "CG", order)).interpolate(mesh.coordinates)
            for mesh in mh]
        ho_meshes = [Mesh(T) for T in Ts]
        from collections import defaultdict
        for i, m in enumerate(ho_meshes):
            m._shared_data_cache = defaultdict(dict)
            for k in mh[i]._shared_data_cache:
                if k != "hierarchy_physical_node_locations":
                    m._shared_data_cache[k] = mh[i]._shared_data_cache[k]
        mh = HierarchyBase(
            ho_meshes, mh.coarse_to_fine_cells,
            mh.fine_to_coarse_cells,
            refinements_per_level=mh.refinements_per_level, nested=mh.nested
        )
        if project_refinements_to_cad:
            for mesh in mh:
                project_to_cad(mesh, cad)
        else:
            project_to_cad(mh[0], cad)
            for i in range(1, len(mh)):
                prolong(Ts[i-1], Ts[i])
    return mh
    
def project_mesh_to_cad_2d(mesh, cad):

    from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
    from OCC.Core.gp import gp_Pnt
    from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnCurve

    coorddata = mesh.coordinates.dat.data
    ids = mesh.exterior_facets.unique_markers

    filt = lambda arr: arr[numpy.where(arr < mesh.coordinates.dof_dset.size)[0]]
    boundary_nodes = {id: filt(mesh.coordinates.function_space().boundary_nodes(int(id), "topological")) for id in ids}

    for id_ in ids:  
        for node in boundary_nodes[id_]:
            #print(node)
            coords = coorddata[node,:]
            best_coords = coords
            dist_old = np.inf
            for edge in cad.edges():
                curve = BRepAdaptor_Curve(edge)
                pt = gp_Pnt(*coorddata[node, :],0)
                proj = GeomAPI_ProjectPointOnCurve(pt, curve.Curve().Curve())
                if proj.NbPoints() > 0:
                    projpt = proj.NearestPoint()
                    projected_coords = np.array(projpt.Coord()[0:2])
                    dist = np.linalg.norm(coords-projected_coords)
                    if dist_old > dist:
                        best_coords = projected_coords
                        dist_old = dist
                    
                    #print(coords, projected_coords, np.linalg.norm(coords-projected_coords))    
            #print(distances)
            coorddata[node,:] = best_coords
    return

def other():
    for (id, edge) in zip(ids, cad.edges()):
        owned_nodes = boundary_nodes[id]
        for other_id in ids:
            if id == other_id:
                continue
            owned_nodes = numpy.setdiff1d(owned_nodes, boundary_nodes[other_id])

        curve = BRepAdaptor_Curve(edge)

        for node in owned_nodes:
            pt = gp_Pnt(*coorddata[node, :], 0)
            proj = GeomAPI_ProjectPointOnCurve(pt, curve.Curve().Curve())
            if proj.NbPoints() > 0:
                projpt = proj.NearestPoint()
                coorddata[node, :] = projpt.Coord()[0:2]
            else:
                warnings.warn("Projection of point %s onto curve failed" % coorddata[node, :])


if __name__=="__main__":
    mesh = Mesh("sphere.msh")
    mh = create_labeled_hierarchy(mesh,"sphere.step", levels = 1 ,order = 2)
    for i,m in enumerate(mh):
        render(m,i)

