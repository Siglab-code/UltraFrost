import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot') 
import transfer_function
import clas

# load experimental transfer function data
da = np.loadtxt('measurement5.csv') 
y_data = da[:,1]/np.max(np.abs(da[:,1]))
f2 =  da[:,0]
loss1 = []
para = []
disp_curve = []
vp_m = 1976
vpp = []

# define aim function 
def aimfuc(x): 
    global f2, y_data, vp_m 
    para.append(x)
    rho1 =np.array([2650]*1) 
    H =np.array([0.17])
    KKs = [x[0]*10**9] 
    muus = [x[1]*10**9] 
    porosity = np.array([x[2]]) 
    Sr =np.array([x[3]])
    phiw1 = porosity*Sr 
    phii1 = porosity*(1-Sr) 
    n1 = 1; layer =2; rmax = 0.1/2; r = 0; th = 0.17 
    # calculate theoretical transfer function and P1 wave velocity
    yt, vp = transfer_function.frozen(f2, KKs, muus, rho1, phiw1, phii1,  n1, layer, r, rmax, th)   
    vpp.append(vp)
    yt = np.array(yt)/np.max(np.abs(yt)) 
    # Monitoring resutls
    # plt.clf() 
    # disp_curve.append(yt)   
    # plotyt = np.array(disp_curve)
    loss = np.abs(yt - y_data)
    los_func = np.abs(vp-vp_m)+ 1000*np.sqrt(np.sum(loss**2))
    # loss1.append(los_func)  
    # loc = np.where(np.array(loss1) == np.array(loss1).min())[0] 
    # plt.clf() 
    # plt.plot(f2, y_data, label= 'measurement')  
    # plt.plot(f2, plotyt[loc[0], :], ':', label='Prediction')   
    # plt.legend()
    # plt.savefig('dispersion_update.png')   
    return los_func

# Run the optimization  
for i in range(9):
    i = i 
    Ei = (6, 10)   # bulk modulus (GPa)
    mui = (5,6)    # shear modulus (GPa)
    po = (0.45,0.55)  # porosity
    sri = (0.05+(i) * 0.1, 0.05+(i + 1) * 0.1)  # degree of saturation of unfrozen water
    print(sri)
    srch = clas.Searcher(
        objective= aimfuc,
        limits=[Ei, mui, po, sri],
        num_samp=20,
        num_resamp=10,
        maximize=False,
        verbose=True
        ) 
    a, b = srch.update(50)

