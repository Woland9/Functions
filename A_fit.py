# Frequency
def a_freq_xy(f,PSDX,PSDY,lims):
    idxf01 = np.argmin(abs(f-lims[0]));
    idxf02 = np.argmin(abs(f-lims[1]));
    idxfs = list(linspace(idxf01,idxf02,idxf02-idxf01+1,dtype=int));
    idxfsok = idxfs;
    mult = array([1,2,3]);
    for idxt in range(len(mult)):
        dfexc = 3000/mult[idxt];
        try:
            idxf01 = np.argmin(abs(f-mult[idxt]*item+dfexc));
            idxf02 = np.argmin(abs(f-mult[idxt]*item-dfexc));
            idxfsexc = list(linspace(idxf01,idxf02,idxf02-idxf01+1,dtype=int));
            idxfsok = setdiff1d(idxfsok,idxfsexc);
        except:
            print('')    
    f1 = f[idxfsok][argmax(PSDX[idxfsok])];
    f2 = f[idxfsok][argmax(PSDY[idxfsok])];
    if(f1-f2<250):
        wtemp = f[idxfsok]*2*pi;
        PSDtemp = PSDX[idxfsok];
        (dump0,par_out,dump1,dump2) = a_fit1(wtemp/2/pi,PSDtemp,f1,1.5e4,g0);
        f2 = wtemp[np.argmax(PSDtemp-a_model1(par_out,wtemp))]/2/pi;
    fx = f1 if f1>f2 else f2;
    fy = f1 if f1<f2 else f2;
    return(fx,fy);
def a_freq_z(f,PSD,lims):
    idxf01 = np.argmin(abs(f-lims[0]));
    idxf02 = np.argmin(abs(f-lims[1]));
    # fz = integrate.trapz(f[idxf01:idxf02]*PSDX[idxf01:idxf02],f[idxf01:idxf02])/integrate.trapz(PSDX[idxf01:idxf02],f[idxf01:idxf02]);
    fz = f[np.argmax(PSD[idxf01:idxf02])+idxf01]
    return(fz)

# Fit 1 peak
def a_model1(v,x):
    return v['a0']*(v['f0']*v['g0'])**2/((v['f0']**2 - x**2)**2+v['g0']**2*x**2) + funBG(x/2/pi);
def a_fcn2min1(params, x, data):
    v = params.valuesdict()
    return a_model1(v,x) - data;
def a_fit1(f,PSD,f0,fcut,g0):
    wfit = 2*pi*f[np.argmin(abs(f-(f0-fcut))):np.argmin(abs(f-(f0+fcut)))];
    PSDfit = PSD[np.argmin(abs(f-(f0-fcut))):np.argmin(abs(f-(f0+fcut)))];
    params = lmfit.Parameters()
    params.add('a0', value=PSD[np.argmin(abs(f-f0))], min=0, max=PSD[np.argmin(abs(f-f0))]*2)
    params.add('f0', value=2*pi*f0, min=2*pi*f0*.99, max=2*pi*f0*1.01)
    params.add('g0', value=g0, min=g0*0.9, max=g0*1.1)    
    minner = lmfit.Minimizer(a_fcn2min1, params, fcn_args=(wfit, PSDfit));
    result = minner.minimize();
    return(result.init_values,result.params.valuesdict(),wfit/2/pi,PSDfit);

# Fit 2 peaks
def a_model2(v,x):
    return v['a0']*(v['f0']*v['g0'])**2/((v['f0']**2 - x**2)**2+v['g0']**2*x**2) + v['a1']*(v['f1']*v['g1'])**2/((v['f1']**2 - x**2)**2+v['g1']**2*x**2) + funBG(x/2/pi);
def a_fcn2min2(params, x, data):
    v = params.valuesdict()
    return a_model2(v,x) - data;
def a_fit2(f,PSD,f0,f1,fcut,g0):
    wfit = 2*pi*f[np.argmin(abs(f-(min(f0,f1)-fcut))):np.argmin(abs(f-(max(f0,f1)+fcut)))];
    PSDfit = PSD[np.argmin(abs(f-(min(f0,f1)-fcut))):np.argmin(abs(f-(max(f0,f1)+fcut)))];
    params = lmfit.Parameters()
    params.add('a0', value=PSD[np.argmin(abs(f-f0))], min=0, max=PSD[np.argmin(abs(f-f0))]*2)
    params.add('f0', value=2*pi*f0, min=2*pi*f0*.99, max=2*pi*f0*1.01)
    params.add('g0', value=g0, min=g0*0.99, max=g0*1.01)
    params.add('a1', value=PSD[np.argmin(abs(f-f1))], min=0, max=PSD[np.argmin(abs(f-f1))]*2)
    params.add('f1', value=2*pi*f1, min=2*pi*f1*.99, max=2*pi*f1*1.01)
    params.add('g1', value=g0, min=g0*0.99, max=g0*1.01)  
    minner = lmfit.Minimizer(a_fcn2min2, params, fcn_args=(wfit, PSDfit));
    result = minner.minimize();
    return(result.init_values,result.params.valuesdict(),wfit/2/pi,PSDfit);

# Fit 2 peaks + Z-leaks
def a_model2z(v,x):
    return (v['a0']*(v['f0']*v['g0'])**2/((v['f0']**2 - x**2)**2+v['g0']**2*x**2) + v['a1']*(v['f1']*v['g1'])**2/((v['f1']**2 - x**2)**2+v['g1']**2*x**2) + 
            v['azm']*(v['fzm']*v['gzm'])**2/((v['fzm']**2 - x**2)**2+v['gzm']**2*x**2) + v['azp']*(v['fzp']*v['gzp'])**2/((v['fzp']**2 - x**2)**2+v['gzp']**2*x**2) + funBG(x/2/pi));
def a_fcn2min2z(params, x, data):
    v = params.valuesdict()
    return a_model2z(v,x) - data;
def a_fit2z(f,PSD,f0,f1,fcut,g0):        
    wfit = 2*pi*f[np.argmin(abs(f-(min(f0,f1)-fcut))):np.argmin(abs(f-(max(f0,f1)+fcut)))];
    PSDfit = PSD[np.argmin(abs(f-(min(f0,f1)-fcut))):np.argmin(abs(f-(max(f0,f1)+fcut)))];
    params = lmfit.Parameters()
    params.add('a0', value=PSD[np.argmin(abs(f-f0))], min=0, max=PSD[np.argmin(abs(f-f0))]*2)
    params.add('f0', value=2*pi*f0, min=2*pi*f0*.99, max=2*pi*f0*1.01)
    params.add('g0', value=g0, min=g0*0.99, max=g0*1.01)
    params.add('a1', value=PSD[np.argmin(abs(f-f1))], min=0, max=PSD[np.argmin(abs(f-f1))]*2)
    params.add('f1', value=2*pi*f1, min=2*pi*f1*.99, max=2*pi*f1*1.01)
    params.add('g1', value=g0, min=g0*0.99, max=g0*1.01)
    fzm = f0 - fz;
    fzp = f0 + fz;
    params.add('azm', value=PSD[np.argmin(abs(f-fzm))], min=0, max=PSD[np.argmin(abs(f-fzm))]*2)
    params.add('fzm', value=2*pi*fzm, min=2*pi*fzm*.99, max=2*pi*fzm*1.01)
    params.add('gzm', value=g0, min=g0/2, max=g0*2)
    params.add('azp', value=PSD[np.argmin(abs(f-fzp))], min=0, max=PSD[np.argmin(abs(f-fzp))]*2)
    params.add('fzp', value=2*pi*fzp, min=2*pi*fzp*.99, max=2*pi*fzp*1.01)
    params.add('gzp', value=g0, min=g0/2, max=g0*2)  
    minner = lmfit.Minimizer(a_fcn2min2z, params, fcn_args=(wfit, PSDfit));
    result = minner.minimize();
    return(result.init_values,result.params.valuesdict(),wfit/2/pi,PSDfit);

# Refitting
def a_refit2_1(f,PSD,par_in,par_out,lims,params):
    par_in2 = {}; par_out2 = {};
    if par_out['X']['a0']/par_out['X']['a1']>1e4:
        print('Refitting')
        (par_in2['X'],par_out2['X'],dump,dump) = a_fit1(f,PSDX,par_in['X']['f0']/2/pi,lims['cut']*par_in['Z']['f0']/2/pi,params['g0']);
        par_in['X']['a1'] = 0; par_out['X']['a1'] = par_in['X']['a1'];
        par_in['X']['f1'] = 2*pi*fx; par_out['X']['f1'] = par_in['X']['f1'];
        par_in['X']['g1'] = params['g0']; par_out['X']['g1'] = par_in['X']['g1'];
    if par_out['X']['a1']/par_out['X']['a0']>1e4:
        print('Refitting')
        (par_in2['X'],par_out2['X'],dump,dump2) = a_fit1(f,PSDX,par_in['X']['f1']/2/pi,lims['cut']*par_in['Z']['f0']/2/pi,params['g0']);
        par_in['X']['a0'] = 0; par_out['X']['a0'] = par_in['X']['a0'];
        par_in['X']['f0'] = 2*pi*fx; par_out['X']['f0'] = par_in['X']['f0'];
        par_in['X']['g0'] = params['g0']; par_out['X']['g0'] = par_in['X']['g0'];
        par_in['X']['a1'] = par_in2['X']['a0']; par_out['X']['a1'] = par_out2['X']['a0'];
        par_in['X']['f1'] = par_in2['X']['f0']; par_out['X']['f1'] = par_out2['X']['f0'];
        par_in['X']['g1'] = par_in2['X']['g0']; par_out['X']['g1'] = par_out2['X']['g0'];
    if par_out['X']['a0']/par_out['Y']['a1']>1e4:
        print('Refitting')
        (par_in2['Y'],par_out2['Y'],dump,dump) = a_fit1(f,PSDY,par_in['Y']['f0']/2/pi,lims['cut']*par_in['Z']['f0']/2/pi,params['g0']);
        par_in['Y']['a1'] = 0; par_out['Y']['a1'] = par_in['Y']['a1'];
        par_in['Y']['f1'] = 2*pi*fx; par_out['Y']['f1'] = par_in['Y']['f1'];
        par_in['Y']['g1'] = params['g0']; par_out['Y']['g1'] = par_in['Y']['g1'];
    if par_out['Y']['a1']/par_out['Y']['a0']>1e4:        
        print('Refitting')
        (par_in2['Y'],par_out2['Y'],dump,dump) = a_fit1(f,PSDY,par_in['Y']['f1']/2/pi,lims['cut']*par_in['Z']['f0']/2/pi,params['g0']);
        par_in['Y']['a0'] = 0; par_out['Y']['a0'] = par_in['Y']['a0'];
        par_in['Y']['f0'] = 2*pi*fx; par_out['Y']['f0'] = par_in['Y']['f0'];
        par_in['Y']['g0'] = params['g0']; par_out['Y']['g0'] = par_in['Y']['g0'];
        par_in['Y']['a1'] = par_in2['Y']['a0']; par_out['Y']['a1'] = par_out2['Y']['a0'];
        par_in['Y']['f1'] = par_in2['Y']['f0']; par_out['Y']['f1'] = par_out2['Y']['f0'];
        par_in['Y']['g1'] = par_in2['Y']['g0']; par_out['Y']['g1'] = par_out2['Y']['g0'];
    return(par_in,par_out)
