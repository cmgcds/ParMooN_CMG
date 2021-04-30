import os
import numpy as np
import pandas as pd
import mpmath as mp

def getPe(b, epsilon, h):
    return (abs(b)*h)/(2. * epsilon);


def getTau(b, epsilon, h):
    Pe = getPe(b,epsilon, h);
    tau = h/(2.*abs(b)) * (mp.coth(Pe) - 1./Pe);
    return tau;

def generateDataset(size, filename):
    bRange = np.linspace(0.0001, 10, size);
    eRange = np.linspace(1e-12, 10, size);
    hRange = np.linspace(1e-6, 1, size);

    data = np.zeros(shape=(size**3,4));
    count = 0;
    for b in bRange:
        for e in eRange:
            for h in hRange:
                data[count,0] = b;
                data[count,1] = e;
                data[count,2] = h;
                data[count,3] = getTau(b,e,h);
                count += 1;
                pass;
            pass;
        pass;

    np.savetxt(filename, data);
    pass;


if __name__ == "__main__":
    generateDataset(10,"data1.dat");
    #print(getTau(1.250602846,0.010784235,0.000462782))
    pass;





    

