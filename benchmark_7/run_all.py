from hyperion.model import Model
from cli            import command_line_arguments, handle_command_line_arguments, filename
from plot_results   import plot_results
from setup_model    import setup_model
from export_to_fits import export_to_fits

cli = command_line_arguments()

# Print comments to screen:
cli.verbose             = 1

# Resolution in pc:
cli.resolution          = 0.1

# Photons per grid cell:
cli.photons_temperature = 1
cli.photons_raytracing  = 1
cli.photons_imaging     = 1

i = 0
for cli.filament in ["linear", "spiraling"]:
	for cli.sources in ["external", "stellar"]:
		for cli.opticaldepth in [5.0, 20.0, 100.0]:
			for cli.mode in ["temperature", "images", "sed"]:
		
				# Setup step:
				setup_model(cli);
		
				# Compute step:
				file = filename(cli, cli.mode)
				model = Model.read(file+".rtin")
				model.write("temp.rtin")
				model.run(filename=file+".rtout", logfile=file+".hyout", overwrite=True)

				# Plotting step:
				#TODO: if(cli.mode == "temperature"): visualize_setup + visualize_results => plot_setup
				if(cli.mode != "temperature"):
					plot_results(cli);
			
				# Export in fits format:
				if(cli.mode == "images"):
					export_to_fits(cli);


			i += 1
			print("Model number " + str(i) + " out of 12 completed.")
