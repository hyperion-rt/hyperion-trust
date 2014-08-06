import sys, getopt


class command_line_arguments(object):
    __slots__ = ['filament', 'sources', 'opticaldepth', 'mode', 'resolution', 'photons_temperature', 'photons_raytracing', 'photons_imaging', 'verbose']


def handle_command_line_arguments(cli):
	Usage_str = \
		"Usage: python " + sys.argv[0] + " <args> \n\
		with mandatory <args>: \n\
		--filament            <'linear'   (1) or 'spiraling' (2)> \n\
		--sources             <'external' (a) or 'stellar'   (b)> \n\
		--opticaldepth        <optical depth of the filament at 1 micron> \n\
		--mode                <'temperature' (t), 'images' (i) or 'sed' (s)> \n\
		--resolution          <resolution in parsec> \n\
		--photons             <number of photons per cell> \n\
		and optional <args>: \n\
		--photons_temperature <number of photons per cell, only used for temperature computation> \n\
		--photons_raytracing  <number of photons per cell, only used for ray-tracing> \n\
		--photons_imaging     <number of photons per cell, only used for imaging and SEDs> \n\
		--verbose             <prints comments on physical, numerical, and imgaging setup>"
	
	#
	# Initialize command-line arguments:
	#
	cli.filament            = ''
	cli.sources             = ''
	cli.opticaldepth        = -1.0
	cli.mode                = ''
	cli.resolution          = -1.0
	photons                 = -1
	cli.photons_temperature = -1
	cli.photons_raytracing  = -1
	cli.photons_imaging     = -1
	cli.verbose             = 0
	
	#
	# Read command-line arguments:
	#
	if(len(sys.argv)) < 13:
		print("ERROR: Too few command-line arguments were given.")
		print(Usage_str)
		sys.exit(2)
	try:
		opts, args = getopt.getopt(sys.argv[1:],"f:s:o:m:r:p:pt:pr:pi:v",["filament=","sources=","opticaldepth=","mode=","resolution=","photons=","photons_temperature=","photons_raytracing=","photons_imaging=","verbose"])
	except getopt.GetoptError:
		print("ERROR: Reading command-line arguments failed.")
		print(Usage_str)
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print(Usage_str)
			sys.exit()
		elif opt in ("-o", "--opticaldepth"):
			cli.opticaldepth = float(arg)
		elif opt in ("-f", "--filament"):
			cli.filament = arg
		elif opt in ("-s", "--sources"):
			cli.sources = arg
		elif opt in ("-m", "--mode"):
			cli.mode = arg
		elif opt in ("-r", "--resolution"):
			cli.resolution = float(arg)
		elif opt in ("-p", "--photons"):
			photons = int(arg)
			cli.photons_temperature = photons
			cli.photons_raytracing  = photons
			cli.photons_imaging     = photons
		elif opt in ("-pt", "--photons_temperature"):
			cli.photons_temperature = int(arg)
		elif opt in ("-pr", "--photons_raytracing"):
			cli.photons_raytracing = int(arg)
		elif opt in ("-pi", "--photons_imaging"):
			cli.photons_imaging = int(arg)
		elif opt in ("-v", "--verbose"):
			cli.verbose = 1
	
	#
	# Abbreviations of input parameters:
	#
	if(cli.filament == "1"): cli.filament = "linear"
	if(cli.filament == "2"): cli.filament = "spiraling"
	if(cli.sources  == "a"): cli.sources  = "external"
	if(cli.sources  == "b"): cli.sources  = "stellar"
	if(cli.mode     == "t"): cli.mode     = "temperature"
	if(cli.mode     == "i"): cli.mode     = "images"
	if(cli.mode     == "s"): cli.mode     = "sed"
	
	#
	# Check command-line arguments:
	#
	if(cli.opticaldepth <= 0):
		print(Usage_str)
		sys.exit(2)
	if(cli.filament != 'linear' and cli.filament != 'spiraling'):
		print(Usage_str)
		sys.exit(2)
	if(cli.sources != 'external' and cli.sources != 'stellar'):
		print(Usage_str)
		sys.exit(2)
	if(cli.mode != 'temperature' and cli.mode != 'images' and cli.mode != 'sed'):
		print(Usage_str)
		sys.exit(2)
	if(cli.resolution <= 0):
		print(Usage_str)
		sys.exit(2)
	if(cli.photons_temperature <= 0):
		print(Usage_str)
		sys.exit(2)
	if(cli.photons_raytracing <= 0):
		print(Usage_str)
		sys.exit(2)
	if(cli.photons_imaging <= 0):
		print(Usage_str)
		sys.exit(2)
	
	#
	# Print command-line arguments:
	#
	if(cli.verbose):
		print("Physical setup:")
		print(" The filament configuration is", cli.filament)
		print(" The radiation source is", cli.sources)
		print(" The optical depth of the filament at 1 micron is", cli.opticaldepth)
	return


def filename(cli, access):
	
	file = "./data/"
	#file = "./test/"
	file += "TRUST-VII" + "_filament=" + cli.filament + "_sources=" + cli.sources + "_opticaldepth=" + str(cli.opticaldepth)
	
	if(access == "temperature"):
		file += "_mode=" + "temperature"
	else:
		file += "_mode=" + cli.mode

	file += "_resolution=" + str(cli.resolution) + "pc" + "_photons_temperature=" + str(cli.photons_temperature) + "percell"

	if(access != "temperature"):
		file += "_photons_raytracing=" + str(cli.photons_raytracing) + "percell" + "_photons_imaging=" + str(cli.photons_imaging) + "percell"
		
	return file
