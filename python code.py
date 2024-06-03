import math
import numpy
numpy.set_printoptions (3, suppress=True)

nodes=int(input('Total number of nodes: '))
elements=int(input('Total number of elements: '))
xc=[]
yc=[]
for i in range (nodes):
    x=float(input('X coordinate of node '+str(i+1)+': '))
    y=float(input('Y coordinate of node '+str(i+1)+': '))
    xc.append(x)
    yc.append(y)

area=float(input("Area of cross section: "))
E=float(input("Modulus of Elasticity: "))
nst=[]
nend=[]
elen=[]
k=[]
ecos=[]
esin=[]
for i in range (elements):
    a=int(input('Start node of element '+str(i+1)+': '))
    b=int(input('End node of element '+str(i+1)+': '))
    k1=float(input('k of '+str(i+1)+':' ))
    x1=float(xc[a-1])
    y1=float(yc[a-1])
    x2 = float(xc[b - 1])
    y2 = float(yc[b - 1])
    l=math.sqrt((x2-x1)**2 + (y2-y1)**2)
    #k1=1/l
    #k1=area*E/l
    cos=(x2-x1)/l
    sin=(y2-y1)/l
    nst.append(a)
    nend.append(b)
    elen.append(l)
    k.append(k1)
    ecos.append(cos)
    esin.append(sin)
#print(elen)
#print(k)
#print(ecos)
#print(esin)

kmat=[]
for i in range(elements):
    cc=float(ecos[i])**2
    ss=float(esin[i])**2
    cs=float(ecos[i])*float(esin[i])
    mat=k[i]*numpy.array([[cc, cs, -cc, -cs],
                    [cs, ss, -cs, -ss],
                    [-cc, -cs, cc, cs],
                    [-cs, -ss, cs, ss]])
    kmat.append(mat)

gmat=[]
for i in range (elements):
    m=nst[i]*2
    n=nend[i]*2
    add=[m-1, m, n-1, n]
    gmat1=numpy.zeros((nodes*2, nodes*2))
    emat=kmat[i]
    for j in range(4):
        for k in range(4):
            a=add[j]-1
            b=add[k]-1
            gmat1[a,b]=emat[j,k]
    gmat.append(gmat1)
gsmat=numpy.zeros((nodes*2, nodes*2))
for mat in gmat:
    gsmat=gsmat+mat
print('\nGlobal Stifness Matrix of the Truss:\n')
print(numpy.around(gsmat, 3))

dl=[]
fl=[]
for i in range (nodes):
    a=str('u')+str(i+1)
    dl.append(a)
    b=str('v')+str(i+1)
    dl.append(b)
    c=str('fx')+str(i+1)
    fl.append(c)
    d=str('fy')+str(i+1)
    fl.append(d)

dis=numpy.ones((nodes*2, 1))
supports=int(input('\nNo of nodes having support: '))
condition=['P = Pinned',
       'H = Horizontal restrained',
       'V = Vertical restrained']
for i in range (supports):
    support=int(input('Node number of support: '))
    for j in condition:
        print(j)
    condi=str(input('Condition of support: '))
    if condi in ['P', 'p']:
        dis[support*2-2, 0]=0
        dis[support*2-1, 0]=0
    elif condi in ['H', 'h']:
        dis[support * 2 - 2, 0] = 0
    elif condi in ['V', 'v']:
        dis[support * 2 - 1, 0] = 0
    else:
        print('Condition Invalid')

force=numpy.zeros((nodes*2, 1))
loaded=int(input('\nNumber of Loaded nodes: '))
for i in range (loaded):
    load= int(input('Node number for loading: '))
    fx=float(input('Fx: '))
    fy=float(input('Fy: '))
    force[load*2-2, 0]=fx
    force[load*2-1, 0]=fy

red=[]
for i in range (nodes*2):
    if dis[i,0]==0:
        red.append(i)
rowred=numpy.delete(gsmat, red, 0)
colred=numpy.delete(rowred, red, 1)
gred=colred
fred=numpy.delete(force, red, 0)
dred=numpy.delete(dis, red, 0)

dresult=numpy.matmul(numpy.linalg.inv(gred), fred)
rin=0
for i in range (nodes*2):
    if dis[i,0]==1:
        dis[i,0]=dresult[rin,0]
        rin=rin+1
fresult=numpy.matmul(gsmat, dis)
print('\nDisplacement Matrix of nodes: ')
print(dis)
print('\nForce Matrix of nodes: ')
print(fresult)

newxco = []
newyco = []
count = 0
for i in range(nodes):
    k = xc[i]+dis[count,0]
    newxco.append(k)
    count = count+1
    l = yc[i]+dis[count,0]
    newyco.append(l)
    count = count+1
newlenofel = []
for i in range(elements):
    a, b = nst[i], nend[i]
    x1 = float(newxco[a-1])
    y1 = float(newyco[a-1])
    x2 = float(newxco[b-1])
    y2 = float(newyco[b-1])
    l = math.sqrt((x2-x1)**2+(y2-y1)**2)
    newlenofel.append(l)
numpy.set_printoptions(3, suppress=False)
elstrain = numpy.zeros((elements,1))
for i in range(elements):
    elstrain[i,0] = (newlenofel[i]-elen[i])/(elen[i])
numpy.set_printoptions(3, suppress=True)
elstress = numpy.zeros((elements, 1))
for i in range(elements):
    elstress[i, 0] = E * elstrain[i, 0]
eforce = numpy.zeros((elements,1))
for i in range(elements):
    eforce[i,0] = area * elstress[i,0]
print('\nForce matrix of elements element')
print(eforce)

