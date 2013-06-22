import numpy

def nonzero(tensor):
    ts = numpy.array(tensor)
    shape_array = ts.shape
    size = [0 for i in range(len(shape_array))]
    for i in range(len(shape_array)):
        size[i] = shape_array[i]
    if len(size) == 4:
        for i in range(size[0]):
            for j in range(size[1]):
                for k in range(size[2]):
                    for l in range(size[3]):
                        if tensor[i][j][k][l] != 0:
                            print "["+str(i)+"]["+str(j)+"]["+str(k)+"]["+str(l)+"] is "+str(tensor[i][j][k][l])+"."
        return
    if len(size) == 3:
        for i in range(size[0]):
            for j in range(size[1]):
                for k in range(size[2]):
                    if tensor[i][j][k] != 0:
                            print "["+str(i)+"]["+str(j)+"]["+str(k)+"] is "+str(tensor[i][j][k])+"."
        return
    if len(size) == 2:
        for i in range(size[0]):
            for j in range(size[1]):
                    if tensor[i][j] != 0:
                            print "["+str(i)+"]["+str(j)+"]is "+str(tensor[i][j])+"."
        return
    else:
        return "error"
                                                                                             
#metric tensor default index structure: _alpha _beta                    
def g(metric):
    components = [0 for i in range(len(metric)^2)]
    i = 0
    while i<len(components):
        components[i + i/len(metric)] = metric[i/len(metric)]
        i += len(metric)
    g = matrix(len(metric),len(metric),components)
    return g

#christoffel symbols tensor default index structure: ^alpha _beta _gamma formula on p. 174 gravity. 
def christoffel(metric):
    N = 4
    gamma = [[[0 for i in range(N)] for j in range(N)] for k in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                gamma[i][j][k] = (((1/2)*( diff(g[i][j],coords[k]) + diff(g[i][k],coords[j]) - diff(g[j][k],coords[i]) ) )*(~g)[i][i]).full_simplify()
    return gamma

#riemann tensor default index structure: ^alpha _beta _gamma _delta. formula p. 452 gravity.
def riemann(metric):
    riemann = [[[[0 for i in range(N)] for j in range(N)] for k in range(N)] for l in range(N)]
    gamma = christoffel(metric)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(N):
                    riemann[i][j][k][l] = diff(gamma[i][j][l],coords[k]) - diff(gamma[i][j][k],coords[l])
                    for o in range(N):
                        riemann[i][j][k][l] += gamma[i][k][o]*gamma[o][j][l] - gamma[i][l][o]*gamma[o][j][k]
                    riemann[i][j][k][l] = (riemann[i][j][k][l]).full_simplify()
    return riemann

#ricci tensor default index structure: _alpha _beta. formula p. 456 gravity.    
def ricci(metric):
    _riemann = riemann(metric)
    ricci = [[0 for i in range(N)] for j in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                ricci[i][j] += _riemann[k][i][k][j]
                ricci[i][j] = (ricci[i][j]).full_simplify()
    return ricci

#einstein tensor default index structure: _alpha _beta. formula p. 483 gravity.
def einstein(metric):
    einstein = [[0 for i in range(N)] for j in range(N)]
    _ricci = ricci(metric)
    R = ricci_scalar(metric)
    for i in range(N):
        for j in range(N):
            einstein[i][j] = (_ricci[i][j] - (1/2)*g[i][j]*R).full_simplify()
    return einstein

def einstein_scalar(metric):
    _einstein = einstein(metric)
    G = 0
    for i in range(N):
        for j in range(N):
            G += ((~g)[i][j]*_einstein[i][j]).full_simplify()
    return G.full_simplify()

def ricci_scalar(metric):
    _ricci = ricci(metric)
    R = 0
    for i in range(N):
        for j in range(N):
            R += (~g)[i][j]*_ricci[i][j]
    return R.full_simplify()

print "GR calculator!"
print "By Rene Zeto, June 2013."

x = var('x')
y = var('y')
z = var('z')
t = var('t')
th = var('th')
phi = var('phi')
rho = var('rho')
p = var('p')
k = var('k')
m = var('m')
f1 = function('f1')
f2 = function('f2')
f3 = function('f3')
f4 = function('f4')

print "Initialized variables:"
print "[x,y,z,t,th,phi,rho,p,k,m]"
print "Initialized functions:"
print "[f1,f2,f3,f4]"
print "This program currently only supports diagonal metrics. Enter functions of multiple values with the values separated by commas (e.g. f1(t-x) as f1(t,x))."
print "When entering your metric components and coordinates, separate the entries with spaces (and nothing else)."
metric_str = raw_input("Enter diagonal metric components: ")
metric = []
for i in range(len(metric_str.split(' '))):
    metric += [SR(metric_str.split(' ')[i])]

coords_str = raw_input("Enter corresponding coordinates, in the same order: ")
coords = []
for i in range(len(coords_str.split(' '))):
    coords += [SR(coords_str.split(' ')[i])]

N = len(metric)
g = g(metric)

print "Your coordinates are:"
print coords
print "Your metric tensor is:" 
print g

print "List of functions (tensor objects are stored as python lists unless otherwise noted):"
print " --nonzero(tensor): lists the nonzero components of the input tensor and the corresponding indicies."
print " --g([diag1, diag2, ...]): takes a list of the diagonal metric components, maps them to a sage matrix for use in calculating the other tensors."
print " --christoffel(g): computes the christoffel symbols, stores them in a list."
print " --riemann(g): computes the riemann curvature tensor."
print " --ricci(g): computes the ricci curvature tensor."
print " --einstein(g): computes the einstein tensor."
print " --ricci_scalar(g): computes the trace of the ricci curvature tensor."
print " --einstein_scalar(g): computes the trace of the einstein tensor."

