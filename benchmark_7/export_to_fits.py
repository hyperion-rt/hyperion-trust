########################################################
## FITS exporter of TRUST-VII filament benchmark test ##
########################################################

#import sys
#import numpy as np
#import math
#import matplotlib.pyplot as plt
from astropy.io import fits
from hyperion.model import ModelOutput
from hyperion.util.constants import pc
from cli import command_line_arguments, handle_command_line_arguments, filename


def export_to_fits(cli):
	
	#
	# Read in the model:
	#
	file = filename(cli, "plot")
	file += ".rtout"
	model = ModelOutput(file)
	
	
	#
	# Write fits file:
	#
	if(cli.mode == "images"):
		
		los = [0 for i in range(3)]
		los[0] = 'x'
		los[1] = 'y'
		los[2] = 'z'

		for k in range(0, 3):
			image = model.get_image(distance=1*pc, units='MJy/sr', inclination=0, component='total', group=k)
			Nwavelength=image.val.shape[2]
			for i in range(0, Nwavelength):
				file = filename(cli, "fits")
				file += "_wavelength=" + str(image.wav[i]) + "micron_los=" + los[k] + ".fits"
				fits.writeto(file, image.val[:, :, i], clobber=True)
				if(cli.verbose):
					print("  The fits file was written to", file)

	else:
		print("ERROR: The specified mode", mode, "is not available. Use 'images' only.")


if(__name__ == "__main__"):
	cli = command_line_arguments()
	handle_command_line_arguments(cli);
	export_to_fits(cli);

