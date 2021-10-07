import numpy as np


def ReadStars(file,Munit=1.65555e+39/2e33,pos=False,recl=4):
    if(recl == 4):
        lint = np.int32
    else:
        assert(recl == 8)
        lint = np.int64
    with open(file,'r') as f:           
        h = np.fromfile(f,dtype=lint,count=1); assert(h[0] == recl)
        n = np.fromfile(f,dtype=lint,count=1)[0]
        h = np.fromfile(f,dtype=lint,count=1); assert(h[0] == recl)
        h = np.fromfile(f,dtype=lint,count=1); assert(h[0] == 24)
        bb = np.fromfile(f,dtype=np.float32,count=6)
        h = np.fromfile(f,dtype=lint,count=1); assert(h[0] == 24)

        p = np.empty((3,n),dtype=np.float64)
        for i in range(3):
            x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 8*n)
            p[i,:] = np.fromfile(f,dtype=np.float64,count=n)
            x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 8*n)
        q = np.empty((4,n),dtype=np.float32)
        for i in range(4):
            x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 4*n)
            q[i,:] = np.fromfile(f,dtype=np.float32,count=n)
            x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 4*n)
        mt = q[0]*Munit
        ms = q[1]*Munit
        ts = q[2]
        zs = q[3]/0.02
        if(pos):
            for i in range(3):
                p[i,:] = (p[i,:]-bb[i])/(bb[i+3]-bb[i])
            return (mt,ms,ts,zs,p)
        else:
            return (mt,ms,ts,zs)


def ReadMesh(file,recl=8):
    if(recl == 4):
        lint = np.int32
    else:
        assert(recl == 8)
        lint = np.int64
    with open(file,'r') as f:
        h = np.fromfile(f,dtype=lint,count=1); assert(h[0] == 12)
        nn = np.fromfile(f,dtype=np.int32,count=3)
        h = np.fromfile(f,dtype=lint,count=1); assert(h[0] == 12)
        n = nn[0]
        assert(n==nn[1] and n==nn[2])

        ret = []
        
        while(1):
            x = np.fromfile(f,dtype=lint,count=1); 
            if(len(x) == 0): return ret
            assert(x[0] == 4*n**3)                
            d = np.fromfile(f,dtype=np.float32,count=n**3).reshape((n,n,n))
            x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 4*n**3)
            ret.append(d)


def ReadProj(file,dx0=6.31947e+22,lev=5,lev_file=6,recl=4):
    if(recl == 4):
        lint = np.int32
        hlen = 5
    else:
        assert(recl == 8)
        lint = np.int64
        hlen = 7
    with open(file,'r') as f:
        h = np.fromfile(f,dtype=np.int32,count=hlen)
        if(recl == 4):
            assert(h[0]==12 and h[4]==12)
            assert(h[1]==h[2] and h[1]==h[2])
            n = h[1]
        else:
            assert(h[0]==12 and h[5]==12)
            assert(h[2]==h[3] and h[2]==h[4])
            n = h[2]
        x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 4*n**3)
        xH2 = np.fromfile(f,dtype=np.float32,count=n**3).reshape((n,n,n))
        x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 4*n**3)
        x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 4*n**3)
        den = np.fromfile(f,dtype=np.float32,count=n**3).reshape((n,n,n))
        x = np.fromfile(f,dtype=lint,count=1); assert(x[0] == 4*n**3)

        w = dx0*0.5**lev_file*1.674e-24/2e33*(3.086e18)**2
        sgas = np.sum(den,axis=0)*w
        smol = np.sum(xH2*den,axis=0)*w
        for l in range(lev,lev_file):
            sgas = 0.25*(sgas[0::2,0::2]+sgas[1::2,0::2]+sgas[0::2,1::2]+sgas[1::2,1::2])
            smol = 0.25*(smol[0::2,0::2]+smol[1::2,0::2]+smol[0::2,1::2]+smol[1::2,1::2])
        return (smol,sgas)
