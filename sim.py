from firedrake import *
import matplotlib.pyplot as plt
from mesh import mesh_hierachy_generator as mhg
import numpy as np
import scipy.io as sio
from subprocess import Popen
import os
import shutil
from intersection import intersection
import scipy.interpolate as si




def mkdir(foldername):
    """
    Making directory at foldername address
    """
    try:
        os.mkdir('./'+foldername)
    except OSError:
        pass
    return

def read_config(lens_config):
    """
    reading and returning a lens_config dictionary. 
    Rr right surface radius of curvature
    Rl left surface radius of curvature
    t thickness of the lens
    h the height of the lens
    bl1 --> extra boundary edge 
    bl2 --> extra boundary edge
    boundaryLayer: if bl1, or bl2 are nonzero, the value of boundryLayer is True
    extruded (True/False) that shows if the lens is extruded (cylindrical lens) or axis-symmetric (spherical lens)
    """
    Rr = lens_config['Rr']
    Rl = lens_config['Rl']
    t = lens_config['t']
    h = lens_config['h']
    bl1 = lens_config['bl1']
    bl2 = lens_config['bl2']
    boundaryLayer = lens_config['boundaryLayer']
    extruded = lens_config['extruded']
    return Rr, Rl, t, h, bl1, bl2, boundaryLayer,extruded

def folderPath(lens_config):
    """ 
    return the folder path corresponding to the lens_config file
    the folder path is as follows
    if concave:
       folder = './output/' + 'Rl' + str(Rl) + '_Rr' + str(Rr) + '_t' + str(t) + '_h' + str(h) + '_blx' + str(bl1) + '_bly' + str(bl2) + '_' + 'cncv_' + meshGeo;
    otherwise
       folder = './output/' + 'Rl' + str(Rl) + '_Rr' + str(Rr) + '_t' + str(t) + '_h' + str(h) + '_blx' + str(bl1) + '_bly' + str(bl2) + '_' + meshGeo;
    """
    
    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config)

    if extruded:
        meshGeo='cyl'
    else:
        meshGeo='sph'

    if boundaryLayer:
        meshGeo += 'BL'
    else:
        bl1=0
        bl2=0

    try:
        concave = lens_config['concave'];
    except KeyError:
        concave = False;

    if concave:
        folder = './output/' + 'Rl' + str(Rl) + '_Rr' + str(Rr) + '_t' + str(t) + '_h' + str(h) + '_blx' + str(
            bl1) + '_bly' + str(bl2) + '_' + 'cncv_' + meshGeo;
    else:
        folder = './output/' + 'Rl' + str(Rl) + '_Rr' + str(Rr) + '_t' + str(t) + '_h' + str(h) + '_blx' + str(
            bl1) + '_bly' + str(bl2) + '_' + meshGeo;

    return folder


def outputPath(lens_config, nu= 0.3,compress = False):
    """ 
    The output path; in the lens trajectory, based on the poisson ratio
    and if it compressing or not
    """
    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config)
    folder = folderPath(lens_config)
    if compress:
        outpath = folder + '/output/cmp_nu_' + str(nu)
    else:
        outpath = folder+'/output/nu_'+str(nu)
    return outpath


def make_mesh(meshGeo, lens_config, meshsize):
    """
    meshGeo is a string either 'extrude' or 'sphere'
                           or  'extrudeBL' 'or sphereBL' 
    lens_config is the confoguration of the lens
    meshsize is the characterestic mesh size
    """
    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config)
    foldername = folderPath(lens_config)
    mkdir(foldername)
    foldername = foldername + '/mesh'
    mkdir(foldername)
    meshFile = foldername + '/mesh.msh';

    try:
        concave = lens_config['concave'];
    except KeyError:
        concave = False;
    if concave:
        if boundaryLayer:
            meshGeoFileName = 'lensBL_con'
        else:
            meshGeoFileName = 'lens_cncv'
    else:
        if boundaryLayer:
            meshGeoFileName = 'lensBL'
        else:
            meshGeoFileName = 'lens'
    # /home/mathiew/Documents/Apps/gmsh-4.6.0-Linux64/bin/
    out = "./gmsh_2 -2 ./mesh/"+meshGeoFileName+".geo -setnumber Rr %f -setnumber Rl %f -setnumber t %f -setnumber h %f -setnumber meshsize %f -setnumber bl1 %f -setnumber bl2 %f -o %s" % (
        Rr, Rl, t, h, meshsize, bl1, bl2, meshFile)
    Popen("./gmsh_2 --version".split())
    gmsh = Popen(out.split())
    gmsh.wait()
    stepFile = foldername + '/'+meshGeo+'_step.geo';
    shutil.copyfile('./mesh/'+meshGeoFileName+'_step.geo', stepFile)
    meshStep = foldername + '/mesh.step'
    out = "./gmsh_2 " + stepFile + " -setnumber Rr %f -setnumber Rl %f -setnumber t %f -setnumber h %f -setnumber bl1 %f -setnumber bl2 %f -parse_and_exit" % (
        Rr, Rl, t, h, bl1, bl2)
    gmsh = Popen(out.split())
    gmsh.wait()
    return meshFile, meshStep

def generate_mesh(lens_config, meshsize):
    """ 
    generates the mesh with right radius Rr, left radius Rl,
    height h , and thickness t, and mesh characterestic length meshsize
    """

    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config)

    mkdir('./output')
    if extruded:
        meshGeo='cyl'
    else:
        meshGeo='sph'

    if boundaryLayer:
        meshGeo += 'BL'
    else:
        # we need to make sure that if the
        # boundarylLayer is false, blx=bly=0
        bl1=0
        bl2=0

    make_mesh(meshGeo, lens_config,meshsize)

    return True


def ggrad(r,x):
    """ 
    calculates the gradient 
    """
    drdR = as_tensor([[r[0].dx(0), 0, r[0].dx(1)],
                      [0, r[0] / x[0], 0],
                      [r[1].dx(0), 0, r[1].dx(1)]])

    return drdR


def defgrad(r,x):
    """
    calculates gradient in the radial coordinates
    """
    drdR = as_tensor([[r[0].dx(0), 0, r[0].dx(1)],
                      [0, r[0] / x[0], 0],
                      [r[1].dx(0), 0, r[1].dx(1)]])
    I = Identity(3)
    F = drdR + I
    return F


def PK1(u,x,extruded = False, E = 1.0, nu=0.3):
    if extruded:
        F = variable(Identity(2)+grad(u))
    else:
        F = variable(defgrad(u,x))
    C = F.T * F
    Ic = tr(C)
    J = det(F)

    mu, lmbda = Constant(E / (2 * (1 + nu))), Constant(E * nu / ((1 + nu) * (1 - 2 * nu)))

    if extruded:
         psi = (mu/2)*(Ic-2)-mu*ln(J) + (lmbda/2)*(ln(J))**2
    else:
        psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J)) ** 2

    return diff(psi, F)

def PK1_inc(u,x,extruded = False, E = 1.0, nu=0.5):
    if extruded:
        F = variable(Identity(2)+grad(u))
    else:
        F = variable(defgrad(u,x))
    C = F.T * F
    Ic = tr(C)
    J = det(F)

    mu = Constant(E / (2 * (1 + nu)))

    if extruded:
         psi = (mu/2)*(Ic-2)
    else:
        psi = (mu / 2) * (Ic - 3) 

    return diff(psi, F)

def sigma_lin(u,E=1,nu =0.499):
    # NOTE this does not work with incompressible

    mu, lmbda = Constant(E / (2 * (1 + nu))), Constant(E * nu / ((1 + nu) * (1 - 2 * nu)))
    eps = sym(grad(u))
    return lmbda*tr(eps)*Identity(2) + 2.0*mu*eps
    

class PointWiseBC(DirichletBC):
    def nodes_to_zero(self, l1, l2):
        import functools
        import numpy as np
        bcnodes = []
        bcnodes1 = []
        b1 = self._function_space.boundary_nodes(l1, self.method)
        b2 = self._function_space.boundary_nodes(l2, self.method)
        bcnodes1.append(b1)
        bcnodes1.append(b2)
        bcnodes1 = functools.reduce(np.intersect1d, bcnodes1)
        bcnodes.append(bcnodes1)
        return bcnodes[0]



def simulate(lens_config, nu=0.3, num_increments = 200, increment_length = 5e-2,compress = False,linear = False):
    """
    the main function that simulates the problem
    """
    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config)
    folder = folderPath(lens_config) + '/mesh/';
    filePathMesh = folder + 'mesh.msh'
    filePathStep = folder + 'mesh.step'
    base = Mesh(filePathMesh)
    # filePathStep = filePathMesh[:-3]+'step'
    mh = mhg.create_labeled_hierarchy(base,filePathStep, levels=0, order=2)
    mesh = mh[-1]

    x = SpatialCoordinate(mesh)
    V = VectorFunctionSpace(mesh, "CG", 2)
    u = Function(V)
    v = TestFunction(V)
    E = 1.0;

    # Storing the strain energy
    Sc = FunctionSpace(mesh, "DG", 1)
    pse = Function(Sc)
    vnms = Function(Sc)
    if extruded:
        Tc = TensorFunctionSpace(mesh, "DG",1)
    else:
        Tc = TensorFunctionSpace(mesh, "DG", 2)
    Green = Function(Tc)
    Green.rename("Strain")
    pse.rename("StrainEnergy")
    vnms.rename("VonMisesStrain")




    # Defining the left and right lens surface
    nleft = V.boundary_nodes(5, 'topological') # left surface of lens marked as 3
    nright = V.boundary_nodes(3, 'topological') # right surface of lens is marked as 5
    coords = interpolate(SpatialCoordinate(mesh), V)
    with coords.dat.vec_ro as xdat:
        pass

    c = xdat.array
    c = c.reshape((-1, 2))
    cl = c[nleft, :]
    cr = c[nright, :]
    # plt.scatter(cr[:, 0], cr[:, 1], c="C0")
    # plt.scatter(cl[:, 0], cl[:, 1], c="C1")
    # plt.show()

    outL = []; # the Left lens surface points
    outR = []; # the right lens surface points

    points_left = cl
    points_right = cr


    newpoints_left = points_left + u.dat.data[nleft, :];
    newpoints_right = points_right + u.dat.data[nright, :];
    outL.append(newpoints_left)
    outR.append(newpoints_right)

    # R.\Omega is the force if the lens is rotating
    if extruded:
        R  = inner(PK1(u,x,extruded = extruded,E=E,nu=nu), grad(v))*dx
    else:
        rhomega = Constant(1e-12)
        R = inner(PK1(u,x,extruded = extruded,E=E,nu=nu), ggrad(v,x)) * x[0] * dx - rhomega * x[0] * v[0] * x[0] * dx

    if linear:
        # This only works for cylindrical case
        R = inner(cauchy_lin(u,E=E),sym(grad(v))*dx

    bc1 = DirichletBC(V.sub(0), Constant(0), 2)
    bc2 = PointWiseBC(V, Constant((0, 0)), 3)

    bc2.nodes = bc2.nodes_to_zero(3,4)
    print(bc2.nodes)
    bcs = [bc1, bc2]

    # b0 = Function(V)
    # b0.interpolate(Constant((0, 1)))
    # nm = VectorSpaceBasis([b0])
    # nm.orthonormalize()

    bcs = [bc1]

    params = {"mat_type": "aij",
              # "snes_monitor": None,
              "snes_rtol": 1e-9,
              "snes_atol": 1e-8,
              "snes_linesearch_type": "l2",
              "ksp_type": "gmres",
              "pc_type": "lu",
              "pc_hypre_type": "boomeramg",
              "pc_factor_mat_solver_type": "mumps"}


    foldername = filePathMesh[:-14] + '/output'
    mkdir(foldername)
    if compress:
        direction = -1;
    else:
        direction = 1;

    if compress:
        foldername = foldername + '/cmp_nu_' + str(nu);
    else:
        foldername = foldername+'/nu_' + str(nu);
    mkdir(foldername)
    outfile = File(foldername + "/result.pvd")
    if linear:
        outfile = File(foldername + "/result_linear.pvd")

    for i in range(num_increments):
        # plt.clf()
        print("Iteration: ", i)

        scale = (i + 1) * increment_length*direction

        bc3 = DirichletBC(V.sub(0), Constant(scale), 4)
        bcs = [bc1, bc2, bc3]

        solve(R == 0, u, bcs, solver_parameters=params)
        newpoints_left = points_left + u.dat.data[nleft, :];
        newpoints_right = points_right + u.dat.data[nright, :];

        outL.append(newpoints_left)
        outR.append(newpoints_right)



        if extruded:
            F = variable(Identity(2)+grad(u))
        else:
            F = variable(defgrad(u,x))
        C = F.T * F
        Ic = tr(C)
        J = det(F)

        mu, lmbda = Constant(E / (2 * (1 + nu))), Constant(E * nu / ((1 + nu) * (1 - 2 * nu)))

        if extruded:
             psi = (mu/2)*(Ic-2)-mu*ln(J) + (lmbda/2)*(ln(J))**2
             Green.interpolate(Constant(0.5) * (C - Identity(2)))
             Strain = Constant(0.5) * (C - Identity(2));
             Straindev = Strain - 1/3.*tr(Strain)*Identity(2);
             SVM = sqrt(2/3.*inner(Straindev,Straindev))
             vnms.interpolate(SVM)
             pse.interpolate(psi)

             outfile.write(u, pse, Green,vnms)
        else:
            psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J)) ** 2
            # Green.interpolate(Constant(0.5) * (C - Identity(3)))
            # pse.interpolate(psi)
            outfile.write(u, pse) # , Green

        # Green.interpolate(Constant(0.5) * (C - Identity(2)))


    np.savez(foldername + '/outLens.npz', outL=outL, outR=outR);
    sio.savemat(foldername + '/outLens.mat',{'outL':outL, 'outR':outR})
    return True


def eig_plus(A):
    return (tr(A) + sqrt(tr(A)**2-4*det(A)))/2


def eig_minus(A):
    return (tr(A) - sqrt(tr(A)**2-4*det(A)))/2 


def eigmatrix(A):
    e1 = eig_plus(A)
    e2 = eig_minus(A)
    
    a = A[0,0]
    b = A[0,1]
    c = A[1,0]
    d = A[1,1]
    tol = 1e-9

    return conditional(ge(abs(c),tol), as_matrix(((e1-d,e2-d), (c,c))),
                       conditional(ge(abs(b), tol),as_matrix(((b,b), (e1-a,e2-a))),
                                   Identity(2)))


def simulate_incompressible(lens_config, nu=0.5, num_increments = 200, increment_length = 5e-2,compress = False):
    """
    the main function that simulates the problem
    """
    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config)
    folder = folderPath(lens_config) + '/mesh/';
    filePathMesh = folder + 'mesh.msh'
    filePathStep = folder + 'mesh.step'
    base = Mesh(filePathMesh)
    # filePathStep = filePathMesh[:-3]+'step'
    mh = mhg.create_labeled_hierarchy(base,filePathStep, levels=0, order=2)
    mesh = mh[-1]

    x = SpatialCoordinate(mesh)
    # Taylor-Hood
    V = VectorFunctionSpace(mesh, "CG", 2)
    Q = FunctionSpace(mesh, "CG", 1)

    W = V*Q
    z = Function(W)
    v = TestFunction(W)
    E = 1.0;

    # Storing the strain energy
    Sc = FunctionSpace(mesh, "DG", 1)
    pse = Function(Sc)
    vnms = Function(Sc)
    ep = Function(Sc, name= "EvP") # Eigenvalues of Strain
    em = Function(Sc, name= "EvM") # Eigenvalues of Strain
    if extruded:
        Tc = TensorFunctionSpace(mesh, "DG",1)
    else:
        Tc = TensorFunctionSpace(mesh, "DG", 2)
    Green = Function(Tc)
    Green.rename("Strain")
    pse.rename("StrainEnergy")
    vnms.rename("VonMisesStrain")
    cauchy = Function(Tc,name="Cauchy")
    vmat = Function(Tc, name="Evecs")




    # Defining the left and right lens surface
    nleft = V.boundary_nodes(5, 'topological') # left surface of lens marked as 3
    nright = V.boundary_nodes(3, 'topological') # right surface of lens is marked as 5
    coords = interpolate(SpatialCoordinate(mesh), V)
    with coords.dat.vec_ro as xdat:
        pass

    c = xdat.array
    c = c.reshape((-1, 2))
    cl = c[nleft, :]
    cr = c[nright, :]
    # plt.scatter(cr[:, 0], cr[:, 1], c="C0")
    # plt.scatter(cl[:, 0], cl[:, 1], c="C1")
    # plt.show()

    outL = []; # the Left lens surface points
    outR = []; # the right lens surface points

    points_left = cl
    points_right = cr

    temp = Function(V)
    newpoints_left = points_left + temp.dat.data[nleft, :];
    newpoints_right = points_right + temp.dat.data[nright, :];
    outL.append(newpoints_left)
    outR.append(newpoints_right)

    # R.\Omega is the force if the lens is rotating
    # extruded=True
    if extruded:
        u,p = split(z)
        F1  = grad(u)+Identity(2)
        J = det(F1)
        Lag_pressure = p*(J-1)*dx
        R1 = derivative(Lag_pressure,z,v) 
        v0,v1 = split(v)
        R  = inner(PK1_inc(u,x,extruded = extruded,E=E), grad(v0))*dx
        R += R1
    else:
        u,p = split(z)
        F1 = defgrad(u,x)
        J = det(F1)
        Lag_pressure = p*(J-1)*x[0]*dx
        R1 = derivative(Lag_pressure,z,v) 
        v0,v1 = split(v)
        R = inner(PK1_inc(u,x,extruded = extruded,E=E,nu=nu), ggrad(v0,x)) * x[0] * dx 
        R += R1

    bc1 = DirichletBC(W.sub(0).sub(0), Constant(0), 2)
    bc2 = PointWiseBC(W.sub(0), Constant((0, 0)), 3)

    bc2.nodes = bc2.nodes_to_zero(3,4)
    print(bc2.nodes)
    bcs = [bc1, bc2]

    bcs = [bc1]

    params = {"mat_type": "aij",
              # "snes_monitor": None,
              "snes_rtol": 1e-9,
              "snes_atol": 1e-8,
              "snes_linesearch_type": "l2",
              "ksp_type": "gmres",
              "pc_type": "lu",
              "pc_factor_mat_solver_type": "mumps"}


    foldername = filePathMesh[:-14] + '/output'
    mkdir(foldername)
    if compress:
        direction = -1;
    else:
        direction = 1;

    if compress:
        foldername = foldername + '/cmp_nu_' + str(nu);
    else:
        foldername = foldername+'/nu_' + str(nu);
    mkdir(foldername)
    outfile = File(foldername + "/result.pvd")

    for i in range(num_increments):
        print("Iteration: ", i)

        scale = (i + 1) * increment_length*direction

        bc3 = DirichletBC(W.sub(0).sub(0), Constant(scale), 4)
        bcs = [bc1, bc2, bc3]

        solve(R == 0, z, bcs, solver_parameters=params)
        u,p = z.split()
        u.rename("Displacement")
        newpoints_left = points_left + u.dat.data[nleft, :];
        newpoints_right = points_right + u.dat.data[nright, :];

        outL.append(newpoints_left)
        outR.append(newpoints_right)



        if extruded:
            F = variable(Identity(2)+grad(u))
        else:
            F = variable(defgrad(u,x))
        C = F.T * F
        Ic = tr(C)
        J = det(F)

        mu = Constant(E / (2 * (1 + nu)))

        if extruded:
             psi = (mu/2)*(Ic-2)
             pse.interpolate(psi)
             Green.interpolate(Constant(0.5) * (C - Identity(2)))
             ep.interpolate(sqrt(eig_plus(C)))   
             em.interpolate(sqrt(eig_minus(C)))   
             vmat.interpolate(eigmatrix(C))

             outfile.write(u, pse, Green,em,ep,vmat)
        else:
            psi = (mu / 2) * (Ic - 3)
            pse.interpolate(psi)
            outfile.write(u, pse) # , Green


    np.savez(foldername + '/outLens.npz', outL=outL, outR=outR);
    sio.savemat(foldername + '/outLens.mat',{'outL':outL, 'outR':outR})
    return True
def sort_points(P,xc):
    P = P[P[:,1].argsort(), 0:2];
    return P;

def focal_aberration_anlaysis(lens_config, n=1.5, nu=0.3, compress=False):

    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config)

    foldername = outputPath(lens_config,nu=nu,compress=compress)
    file = np.load(foldername+'/outLens.npz')

    outL = file['outL']
    outR = file['outR']

    assert len(outL) == len(outR)
    Rot = np.array([[0, -1], [1, 0]])
    smoothing = 0.0;
    distance_right = 200;

    raysX = []; # The xpoints of the rays --> list of arrays --> each array's  row is one ray
    raysY = []; # The ypoints of the rays --> list of arrays --> each array's row is one ray
    focalPoints = []; # the focal points of the rays --> list of arrays
    N_increments = len(outL); # Number of increments
    # raysInitPosY = np.arange(0.01*h, 0.99*h,h/100); # the initial y of the rays
    raysInitPosY = np.arange(0.01*h, 0.8*h,h/100); # the initial y of the rays
    raysInitPosX = -30*np.ones(raysInitPosY.shape); # the initial x of the rays
    raysInitAngle = np.zeros(raysInitPosY.shape); # the initial angle of the rays

    for i in range(N_increments):
        print('Post-processing, {}/{}'.format(i+1,N_increments));


        pointsX = [];
        pointsY = [];

        Pl = outL[i]
        Pr = outR[i]

        Pl = np.matmul(Rot,Pl.T).T;
        Pr = np.matmul(Rot, Pr.T).T;

        xc = (max(Pl[:,0]) + min(Pr[:,0]))/2.0;
        Pl = sort_points(Pl,xc);
        Pr = sort_points(Pr,xc);
        p = np.concatenate((Pl,Pr),axis=0)
        # Pr = np.flipud(Pr);

        Pll = np.insert(Pl, 0, [Pl[0, 0], 2 * Pl[0, 1] - Pl[1, 1]], axis=0)
        Prr = np.insert(Pr, 0, [Pr[0, 0], 2 * Pr[0, 1] - Pr[1, 1]], axis=0)

        left = si.splrep(Pl[:, 1], Pl[:, 0], s=smoothing);
        ynewL = np.arange(0, h, 0.01)
        xnewL = si.splev(ynewL, left, der=0)
        left_der = si.splder(left, n=1)  # spline represenation of the first derivative


        # Fitting a spline to the right lens surface
        right = si.splrep(Pr[:, 1], Pr[:, 0], s=smoothing);
        ynewR = np.arange(0, h, 0.01)
        xnewR = si.splev(ynewR, right, der=0)
        right_der = si.splder(right, n=1)  # spline represantation of the first derivative

        focals = np.zeros(raysInitPosY.shape);

        for ii, yy in enumerate(raysInitPosY):

            p1 = np.array([raysInitPosX[ii], raysInitPosY[ii]]);
            th1 = raysInitAngle[ii];

            p2 = np.array([max(Pl[:,0]), p1[1] + (max(Pl[:,0]) - p1[0]) * np.tan(th1)]);
            xL, yL = intersection(xnewL, ynewL, np.array([p1[0], p2[0]]), np.array([p1[1], p2[1]]))
            normalL = np.array([si.splev(yL[0], left_der), 1]);
            normalL = Rot.dot(normalL);
            normalL = normalL / np.linalg.norm(normalL);

            dir1 = np.array([normalL[0], normalL[1], 0])
            dir2 = np.array([p1[0] - xL[0], p1[1] - yL[0], 0])
            directL = np.cross(dir1, dir2)

            thinL = np.pi - np.arccos(np.dot(normalL, np.array([np.cos(th1), np.sin(th1)])));
            thoutL = np.arcsin(np.sin(thinL) / n);
            ThNegNormL = np.arctan2(-normalL[1], -normalL[0]);
            # thoutLh = ThNegNormL - np.sign(ThNegNormL) * thoutL;
            thoutLh = ThNegNormL + np.sign(directL[2])  * thoutL;

            xR, yR = intersection(xnewR, ynewR, np.array([xL[0], xL[0] + 3 * t]),
                                  np.array([yL[0], yL[0] + 3 * t * np.tan(thoutLh)]))
            normalR = np.array([si.splev(yR[0], right_der), 1]);
            normalR = Rot.dot(normalR);
            normalR = normalR / np.linalg.norm(normalR);

            dir1 = np.array([normalR[0], normalR[1], 0])
            dir2 = np.array([xL[0] - xR[0], yL[0] - yR[0], 0])
            directR = np.cross(dir1, dir2)

            thinR = np.pi - np.arccos(np.dot(normalR, np.array([np.cos(thoutLh), np.sin(thoutLh)])));
            thoutR = np.arcsin(n * np.sin(thinR));
            ThNegNormR = np.arctan2(-normalR[1], -normalR[0]);

            # thoutRh = ThNegNormR - np.sign(ThNegNormR) * thoutR;
            thoutRh = ThNegNormR + np.sign(directR[2]) * thoutR;

            focals[ii] = xR[0] - yR[0] / np.tan(thoutRh);
            pointsX.append(np.array([p1[0], xL[0], xR[0], xR[0]+ distance_right]));
            pointsY.append(np.array([p1[1], yL[0], yR[0], yR[0] + distance_right * np.tan(thoutRh)]));


        raysX.append(np.array(pointsX))
        raysY.append(np.array(pointsY))
        focalPoints.append(focals)

    np.savez(foldername + '/post_process.npz', raysX=raysX, raysY=raysY, focalPoints=focalPoints,
             N_increments=N_increments, raysInitPosY=raysInitPosY, raysInitPosX=raysInitPosX,
             raysInitAngle=raysInitAngle);

    sio.savemat(foldername + '/post_process.mat', {'raysX':raysX, 'raysY':raysY, 'focalPoints':focalPoints,
             'N_increments':N_increments, 'raysInitPosY':raysInitPosY, 'raysInitPosX':raysInitPosX,
             'raysInitAngle':raysInitAngle});


def lens_simulation(lens_config, meshsize, nu, n, num_increments=200, increment_length = 5e-2,compress = False):

    Rr, Rl, t, h, bl1, bl2, boundaryLayer, extruded = read_config(lens_config);
    if type(nu)==list:
        for nu_i in nu:
            generate_mesh(lens_config, meshsize)
            simulate(lens_config, nu=nu_i, num_increments=num_increments, increment_length=increment_length,compress = compress)
            focal_aberration_anlaysis(lens_config, n=n, nu=nu_i,compress = compress)

    elif type(Rr)==list:
        for Rr_i in Rr:
            lens_config['Rr'] = Rr_i
            generate_mesh(lens_config, meshsize)
            simulate(lens_config, nu=nu, num_increments=num_increments, increment_length=increment_length,compress =  compress)
            focal_aberration_anlaysis(lens_config, n=n, nu=nu,compress = compress)
    elif type(Rl) == list:
        for Rl_i in Rl:
            lens_config['Rl'] = Rl_i
            generate_mesh(lens_config, meshsize)
            simulate(lens_config, nu=nu, num_increments=num_increments, increment_length=increment_length,compress = compress)
            focal_aberration_anlaysis(lens_config, n=n, nu=nu,compress = compress)
    elif type(h)==list:
        for h_i in h:
            lens_config['h'] = h_i
            generate_mesh(lens_config, meshsize)
            simulate(lens_config, nu=nu, num_increments=num_increments, increment_length=increment_length,compress = compress)
            focal_aberration_anlaysis(lens_config, n=n, nu=nu,compress = compress)
    elif type(t) == list:
        for t_i in t:
            lens_config['t'] = t_i
            generate_mesh(lens_config, meshsize)
            simulate(lens_config, nu=nu, num_increments=num_increments, increment_length=increment_length,compress = compress)
            focal_aberration_anlaysis(lens_config, n=n, nu=nu,compress = compress)
    elif type(bl1) == list and type(bl2)==list:
        for bl1_i,bl2_i in zip(bl1,bl2):
            lens_config['bl1']= bl1_i
            lens_config['bl2']= bl2_i
            generate_mesh(lens_config, meshsize)
            simulate(lens_config, nu=nu, num_increments=num_increments, increment_length=increment_length,compress = compress)
            focal_aberration_anlaysis(lens_config, n=n, nu=nu,compress = compress)
    else:
        print("Generating Incompressible example")
        generate_mesh(lens_config, meshsize)
        print("DONE")
        simulate_incompressible(lens_config, nu=nu, num_increments=num_increments, increment_length=increment_length,compress = compress)
        focal_aberration_anlaysis(lens_config, n=n, nu=nu,compress = compress)



def main():
    
    meshsize = 0.5;
    n = 1.4;
    nu = 0.3
    increment_length = 5e-2
    num_increments = 200
    compress = False

    lens_config = {'Rr':40, 'Rl':1500, 't':20, 'h':25, 'bl1': 0, 'bl2':0, 'boundaryLayer':False, 'extruded':True,'concave':False}
    #num_increments = 50
    #lens_simulation(lens_config, meshsize, [0.31,0.451], n, num_increments = num_increments, increment_length = increment_length, compress = True)
    num_increments = 100
    lens_simulation(lens_config, meshsize, 0.5, n, num_increments = num_increments, increment_length = increment_length, compress = False)



    #lens_config = {'Rr':40, 'Rl':1500, 't':18, 'h':25, 'bl1': 0, 'bl2':0, 'boundaryLayer':False, 'extruded':False, 'concave':False}
    #num_increments = 50
    #lens_simulation(lens_config, meshsize, [0.31,0.451], n, num_increments = num_increments, increment_length = increment_length, compress = True)
    #num_increments = 100
    #lens_simulation(lens_config, meshsize, [0.31,0.451], n, num_increments = num_increments, increment_length = increment_length, compress = False)


 

if __name__=="__main__":

    main()


