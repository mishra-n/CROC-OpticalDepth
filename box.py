from universe import *
import ifrit as ifrit

class Boxes:
    def __init__(self, Universe, path, resolution=1024, box_size=40):
        
        self.Universe = Universe
        self.path = path
        self.resolution = resolution #pixels
        self.box_size = box_size #cMpc
        self.pixel_size = self.box_size/self.resolution*3.086e24
        
        #list files
        self.file_list = []
        for file in sorted(os.listdir(path)):
            if file.startswith("xd.a=0."):
                self.file_list.append(file)
        
        #loop through directory and list the scale factors
        self.a = np.zeros(len(self.file_list))
        for i, file in enumerate(self.file_list):
            self.a[i] = file[5:11]
         
        self.z = 1/self.a - 1    
            
        #create angular size conversion function (from     
        scales = np.array([0.015, 0.0253, 0.0450, 0.0622, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18])
        angular_size = np.array([32.15, 33.83, 36.41, 38.47, 40.54, 42.83, 45.13, 47.48, 49.9, 52.4])
        
        self.angular = interp.interp1d(scales, angular_size)
        
        
    def angular_size_cut(self, a, closest_a, twoD_sheet):
    
        angular_close = self.angular(closest_a)
        angular_far = self.angular(a)

        pixel_cut = angular_close/ angular_far * self.resolution#twoD_sheet.shape[0]

        subsheet = twoD_sheet[0:int(np.round(pixel_cut)), 0:int(np.round(pixel_cut))]

        x = np.linspace(0,self.resolution, int(np.round(pixel_cut)))
        y = np.linspace(0,self.resolution, int(np.round(pixel_cut)))

        xx, yy = np.meshgrid(x, y)

        f = interp.interp2d(x, y, subsheet)

        xnew = np.arange(self.resolution)
        ynew = np.arange(self.resolution)

        newsheet = f(xnew,ynew)
        
        return newsheet
    
    def bin_box(self, box, binning=1):
        
        if binning==1:
            return box
        else:
            shape = box.shape
            new_shape = shape/binning
            
            box = box.reshape(new_shape[0], shape[0]//new_shape[0], new_shape[0], shape[0]//new_shape[0], new_shape[0], shape[0]//new_shape[0])
            box = box.mean(axis=1).box(axis=2).box(axis=3)

            return box
    
    
    def createUniverse(self, z_max=15, binning=1):
        """
        Unpacks all the data and puts it into numpy arrays [redshift, dim1, dim2, dim3]
        z_max: maximum redshift, saves memory. We don't need data from high redshifts, is little reionizatoin
        binning: must be a multiple of 2
        """
        indices = self.z < z_max
        indexes = np.where(self.z < z_max)[0].tolist()
        
        self.file_list = [self.file_list[index] for index in indexes]
        self.slice_count = sum(indices)
        
        self.z = self.z[indices]
        self.a = self.a[indices]
       
    
        self.n_frac = np.zeros(shape=(len(self.a), self.resolution, self.resolution, self.resolution))
        self.density = np.zeros(shape=(len(self.a), self.resolution, self.resolution, self.resolution))
        
        for i, file in enumerate(self.file_list):
            print(self.a[i], self.z[i])
            self.n_frac[i,:,:,:], self.density[i,:,:,:] = ifrit.ReadMesh(self.path + file)
            
        return self.n_frac, self.density
    
    def binUniverse(self, binning):
        new_res = int(resolution/binning)
        binned_n_frac = np.zeros(shape=(len(self.a), new_res, new_res, new_res))
        binned_density = np.zeros(shape=(len(self.a), new_res, new_res, new_res))
        
        for i in range(n_frac.shape[0]):
            binned_n_frac = self.bin_box(n_frac[i,:,:,:], binning)
            binned_density = self.bin_box(n_density[i,:,:,:], binning)
            
        
        self.n_frac = binned_n_frac
        self.density = binned_density
        self.resolution = new_res
        
        return self.n_frac, self.density
        
    def randOrientUniverse(self):
        np.random.seed(2021)
        
        for i in range(self.slice_count):
            randoms = np.random.randint(0,3,8)
            random_shift = np.random.randint(0,self.resolution,2)

            self.n_frac[i,:,:,:] = np.swapaxes(self.n_frac[i,:,:,:], randoms[0], randoms[1])
            self.density[i,:,:,:] = np.swapaxes(self.density[i,:,:,:], randoms[0], randoms[1])
          #  temp = np.swapaxes(temp, randoms[0], randoms[1])

            self.n_frac = np.flip(self.n_frac[i,:,:,:], randoms[2])
            self.density[i,:,:,:] = np.flip(self.density[i,:,:,:], randoms[2])
           # temp = np.flip(temp, randoms[2])

            self.n_frac[i,:,:,:] = np.roll(self.n_frac[i,:,:,:], shift=random_shift[0], axis=randoms[3])
            self.density[i,:,:,:] = np.roll(self.density[i,:,:,:], shift=random_shift[0], axis=randoms[3])
            
    def tauSlices(self):
        
        self.tau = np.zeros(shape=(len(self.a), resolution, resolution))
        self.dtau = np.zeros(shape=(len(self.a), resolution, resolution))
        self.dtau_patchy = np.zeros(shape=(len(self.a), resolution, resolution))
        self.dtau_fluc = np.zeros(shape=(len(self.a), resolution, resolution))
        
        for i in range(self.slice_count):
            n_b = 1.123e-5 * self.Universe.ombh2 / self.a[i]**3
        self.tau[i,:,:] = self.Universe.sig_T * np.sum((1-self.n_frac[i,:,:,:])*n_b*self.density[i,:,:,:]*a[i]*self.pixel_size, axis=0)            
        self.dtau[i,:,:] = self.Universe.sig_T * np.sum(((1-self.n_frac[i,:,:,:])*self.density[i,:,:,:] - (1-np.mean(self.n_frac[i,:,:,:]))*np.mean(self.density[i,:,:,:]))*n_b*a[i]*self.pixel_size, axis=0)
        
        self.dtau_patchy[i,:,:] = self.Universe.sig_T * np.sum((((1-self.n_frac[i,:,:,:]) - (1-np.mean(self.n_frac[i,:,:,:])))*self.density[i,:,:,:])*n_b*a[i]*self.pixel_size, axis=0)
        
        self.dtau_fluc[i,:,:] = self.Universe.sig_T * np.sum(((1-np.mean(self.n_frac[i,:,:,:]))*(self.density[i,:,:,:] - np.mean(self.density[i,:,:,:])))*n_b*a[i]*self.pixel_size, axis=0)
        
        return self.tau, self.dtau, self.dtau_patchy, self.dtau_fluc
    
    def angularSizeCorrectionTau(self):
        for i in range(self.slice_count):
            self.tau[i,:,:] = angular_size_cut(self.tau[i,:,:], max(self.a), self.a[i])
            self.dtau[i,:,:] = angular_size_cut(self.dtau[i,:,:], max(self.a), self.a[i])
            self.dtau_patchy[i,:,:] = angular_size_cut(self.dtau_patchy[i,:,:], max(self.a), self.a[i])
            self.dtau_fluc[i,:,:] = angular_size_cut(self.dtau_fluc[i,:,:], max(self.a), self.a[i])
        
        return self.tau, self.dtau, self.dtau_patchy, self.dtau_fluc
    
    def genTauMap(self):
        
        self.tau_map = np.sum(self.tau, axis=0)
        self.dtau_patchy_map = np.sum(self.dtau_patchy, axis=0)
        self.dtau_fluc_map = np.sum(self.dtau_fluc, axis=0)
        
        return self.tau_map, self.dtau_patchy_map, self.dtau_fluc_map