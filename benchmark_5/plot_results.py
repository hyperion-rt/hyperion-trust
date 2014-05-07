########################################################
## Visualization of TRUST-VII filament benchmark test ##
########################################################

import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from hyperion.util.constants import pc
from cli import command_line_arguments, handle_command_line_arguments, filename


def plot_results(cli):
    
    file = filename(cli, "plot")
    file += ".rtout"
    
    #
    # Read in the model:
    #
    model = ModelOutput(file)
    
    los = [0 for k in range(3)]
    los[0] = '30degree'
    los[1] = '80degree'
    los[2] = '88degree'
    
    
    if(cli.mode == "images"):
        
        #
        # Extract the quantities
        #
        g = model.get_quantities()
        
        #
        # Get the wall positions:
        #
        ww = g.w_wall / pc
        zw = g.z_wall / pc
        pw = g.p_wall
        
        grid_Nw = len(ww) - 1
        grid_Nz = len(zw) - 1
        grid_Np = len(pw) - 1
        
        #
        # Graphics:
        #
        fig = plt.figure()
        
        Imaxp = [0 for i in range(5)]
        Imaxp[0] = 1e-15 # in W/cm^2
        Imaxp[1] = 1e-14 # in W/cm^2
        Imaxp[2] = 1e-15 # in W/cm^2
        Imaxp[3] = 1e-15 # in W/cm^2
        Imaxp[4] = 1e-18 # in W/cm^2
        
        for k in range(0, 3):
            if(cli.verbose):
                print("Group: ", k)
            
            image = model.get_image(distance=1e+7*pc, units='ergs/cm^2/s', inclination=0, component='total', group=k)
            #source_emit = model.get_image(distance=1e+7*pc, units='MJy/sr', inclination=0, component='source_emit', group=k)
            #dust_emit   = model.get_image(distance=1e+7*pc, units='MJy/sr', inclination=0, component='dust_emit'  , group=k)
            #source_scat = model.get_image(distance=1e+7*pc, units='MJy/sr', inclination=0, component='source_scat', group=k)
            #dust_scat   = model.get_image(distance=1e+7*pc, units='MJy/sr', inclination=0, component='dust_scat'  , group=k)
            
            if(cli.verbose):
                print(" Data cube: ", image.val.shape)
                print(" Wavelengths =", image.wav)
                print(" Uncertainties =", image.unc)
            
            image_Nx=image.val.shape[0]
            image_Ny=image.val.shape[1]
            Nwavelength=image.val.shape[2]
            
            if(cli.verbose):
                print(" Image Nx =", image_Nx)
                print(" Image Ny =", image_Ny)
                print(" Nwavelength =", Nwavelength)
            
            for i in range(0, Nwavelength):
                
                if(cli.verbose):
                    print(" Image #", i,":")
                    print("  Wavelength =", image.wav[i])
                
                image.val[:, :, i] *= 1e-4 # in W/m^2
                
                #Imin = np.min(image.val[:, :, i])
                #Imax = np.max(image.val[:, :, i])
                #Imax = Imaxp[i]
                #Imin = Imax/1e+20
                Imax = np.max(image.val[:, :, i])/5
                Imin = 0.0
                
                if(cli.verbose):
                    print("  Intensity min data values =", np.min(image.val[:, :, i]))
                    print("  Intensity max data values =", np.max(image.val[:, :, i]))
                    print("  Intensity min color-table =", Imin)
                    print("  Intensity max color-table =", Imax)
                
                #ax = fig.add_subplot(2, 1, 2)
                ax = fig.add_subplot(1, 1, 1)
                # 'hot', see http://wiki.scipy.org/Cookbook/Matplotlib/Show_colormaps
                ax.imshow(image.val[:, :, i], vmin=Imin, vmax=Imax, cmap=plt.cm.hot, origin='lower')
                ax.set_xticks([0,100,200,300,400,500], minor=False)
                ax.set_yticks([0,100,200,300,400,500], minor=False)
                ax.set_xlabel('x (pixel)')
                ax.set_ylabel('y (pixel)')
                ax.set_title(str(image.wav[i]) + ' microns' + '\n' + los[k], y=0.88, x=0.5, color='white')
                
                #ax = fig.add_subplot(2, 1, 1)
                #ax.imshow([np.logspace(np.log10(Imin+1e-10),np.log10(Imax/10),100),np.logspace(np.log10(Imin+1e-10),np.log10(Imax/10),100)], vmin=Imin, vmax=Imax/10, cmap=plt.cm.gist_heat)
                #ax.set_xticks(np.logspace(np.log10(Imin+1e-10),np.log10(Imax/10),1), minor=False)
                ##ax.set_xticks(np.linspace(np.log10(Imin+1e-10),np.log10(Imax/10),10), minor=False)
                #ax.set_yticks([], minor=False)
                #ax.set_xlabel('flux (MJy/sr)')
                
                #x = plt.colorbar()
                #print(x)
                
                file = filename(cli, "plot")
                file += "_wavelength=" + str(image.wav[i]) + "micron_los=" + los[k] + ".png"
                
                fig.savefig(file, bbox_inches='tight')
                if(cli.verbose):
                    print("  The image graphics was written to", file)
                plt.clf()
    
    elif(cli.mode == "seds"):
        
        #
        # Graphics:
        #
        fig = plt.figure()
        
        for k in range(0, 3):
            if(cli.verbose):
                print("Group: ", k)
            
            sed = model.get_sed(distance=1e+7*pc, inclination=0, aperture=-1, group=k)
            #units='ergs/cm^2/s' # = default, if distance is specified
            ax = fig.add_subplot(1, 1, 1)
            ax.loglog(sed.wav, sed.val)
            ax.set_xlabel(r'$\lambda$ [$\mu$m]')
            ax.set_ylabel(r'$\lambda F_\lambda$ [ergs/s/cm$^2$]')
            ax.set_xlim(0.09, 1000.0)
            ax.set_ylim(1.e-13, 1.e-7)
            
            file = filename(cli, "plot")
            file += "_los=" + los[k] + ".png"
            fig.savefig(file)
            if(cli.verbose):
                print(" The sed graphics was written to", file)
            plt.clf()
        
        #
        # Data files:
        #
        for k in range(0, 3):
            sed = model.get_sed(distance=1e+7*pc, inclination=0, aperture=-1, group=k)
            file = filename(cli, "plot")
            file += "_los=" + los[k] + ".dat"
            sedtable = open(file, 'w')
            sedtable.write("# wavelength [micron] - flux [erg cm^-2 s^-1]\n")
            for lp in range(0, len(sed.wav)):
                l = len(sed.wav)-lp-1
                line = str("%.4e" % sed.wav[l]) + " " + str("%.4e" % sed.val[l]) + "\n"
                sedtable.write(line)
            sedtable.close()

    else:
        print("ERROR: The specified mode", mode, "is not available. Use 'images' or 'seds' only.")


if(__name__ == "__main__"):
    cli = command_line_arguments()
    handle_command_line_arguments(cli);
    plot_results(cli);

