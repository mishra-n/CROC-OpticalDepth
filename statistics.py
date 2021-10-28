from map import *

class Statistics:
    def __init__(self, Map):
        self.Map = Map
        self.Universe = Map.Universe
        
        self.npix = self.Map.resolution
        self.kfreq = np.fft.fftfreq(self.npix) * self.npix
        self.kfreq2D = np.meshgrid(self.kfreq, self.kfreq)
        self.knrm = np.sqrt(self.kfreq2D[0]**2 + self.kfreq2D[1]**2)
        self.knrm = self.knrm.flatten()
        
        self.kbins = np.arange(0.5, self.npix//2+1, 1.)
        
        self.kvals = 0.5 * (self.kbins[1:] + self.kbins[:-1])
        self.l_1 = 291.9
        self.ell_vals = self.l_1 * self.kvals
        self.lmax = max(self.ell_vals)
        
        
    def crossPowerSpectrum(self, map1, map2):
        map1_k = np.fft.fftn(map1, norm='ortho')
        map2_k = np.fft.fftn(map2, norm='ortho')
        
        self.map1_amplitudes = np.real(map1_k * np.conjugate(map2_k))
        
        Abins, _, _ = stats.binned_statistic(self.knrm, self.map1_amplitudes.flatten(),
                                     statistic = "mean",
                                     bins = self.kbins)
        Abins *= np.pi * (self.kbins[1:]**2 - self.kbins[:-1]**2)
        
        return self.ell_vals, Abins
    
    def autoPowerSpectrum(self, map1):
        
        ell_vals, Abins = self.crossPowerSpectrum(map1, map1)
        
        return self.ell_vals, Abins
    
    
    def theorySpectrum(self, z_range=0):
        
        k, z, P_zk = self.Universe.getCAMBPowerSpectrum(z_range = z_range)
        
        print(k.shape)
        print(z.shape)
        print(P_zk.shape)
        
        comoving_radial_distance = np.vectorize(self.Universe.CAMB_results.comoving_radial_distance)
        
        chi = comoving_radial_distance(z)
        a = 1/(1+z)
        ells = np.logspace(-2, np.log10(self.lmax), 6000)
        
        np.zeros(shape=(len(chi), len(ells)))
        
        file_path = '/data/gnedin/REI/D/Cai.B40.N256L2.sf=1_uv=0.15_bw=10_res=100.WC1/A/rei.log'
        data = np.loadtxt(file_path)
        temp_a = data[:,0]
        temp_z = 1/(temp_a) - 1
        temp_chi = comoving_radial_distance(temp_z)
        XHI_mass = data[:,8]
        XHI_vol = data[:,13]
        print(temp_chi.shape)
        print(XHI_mass.shape)
        
        XHI_vol_func = interp.interp1d(temp_chi, XHI_mass, fill_value=(XHI_mass[0], XHI_mass[-1]), bounds_error=False)
        
        C_ell_integrand=np.zeros(shape=(len(chi), len(ells)))
        
        for d, chi_val in enumerate(chi):
            k_temp = (ells + 1/2)/chi_val
            if chi_val==0:
                continue
            
            Pmm_interp = interp.interp1d(k, P_zk[d,:], fill_value=P_zk[d,-1], bounds_error=False)
            
            C_ell_integrand[d,:] = Pmm_interp(k_temp) / (a[d]**4 * chi_val**2) * XHI_vol_func(chi_val)**2

        
        n_b0 = (0.76+2*0.06) * 1.123e-5 * self.Universe.ombh2 #1/cm**3
        n_b0 = n_b0 * (3.96e24)**2 #1/Mpc**3
        
        sig_T_converted = self.Universe.sig_T / (3.96e24) 
            
        C_ell = n_b0**2 * sig_T_converted**2 * integrate.simps(C_ell_integrand[1::,:], x=chi[1::], axis=0)
        
        return ells, C_ell
            
            
