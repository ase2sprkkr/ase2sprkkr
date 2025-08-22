def main():
    from ase2sprkkr.sprkkr.calculator import SPRKKR
    import os
    import sys
    
    print("Starting SPRKKR calculation...")
    
    try:
        # Get the directory where this script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Set the working directory to the script's directory
        os.chdir(script_dir)
        print(f"Working directory set to: {os.getcwd()}")
        
        # Initialize calculator
        calculator = SPRKKR()
        
        # Set input parameters
        print("Setting up input parameters...")
        calculator.input_parameters = 'torque'
        calculator.input_parameters.CONTROL.DATASET = 'Fe'
        calculator.input_parameters.MODE.MDIR = [1.0, 0.0, 0.0]  # Using floats instead of integers
        calculator.input_parameters.MODE.MALF = 0.0
        calculator.input_parameters.MODE.MBET = 45.0
        calculator.input_parameters.MODE.MGAM = 0.0
        
        # Define potential file path relative to the script location
        potential_file = os.path.abspath('Fe.pot_new')
        
        # Check if potential file exists
        if not os.path.exists(potential_file):
            print(f"\nError: Potential file not found at {potential_file}")
            print("Please ensure the potential file 'Fe.pot_new' is in the same directory as this script.")
            return 1
            
        print(f"Using potential file: {potential_file}")
        print(f"File exists: {os.path.exists(potential_file)}")
        print(f"File size: {os.path.getsize(potential_file) / 1024:.2f} KB")
        
        # Run calculation with explicit path to kkrgen and save output to file
        print("Starting calculation...")
        output_file = os.path.abspath('torque_calculation.out')
        
        try:
            result = calculator.calculate(
                potential=potential_file,
                directory=script_dir,  # Set the working directory for the calculation
                executable_dir='/home/ridha/bin/kkrgen9.4merge',
                print_output=True,
                output_file=output_file
            )
            print(f"Calculation completed! Output saved to {output_file}")
            
            # Read and print the last 20 lines of the output file
            if os.path.exists(output_file):
                print("\n=== Last 20 lines of output ===")
                with open(output_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines[-20:]:
                        print(line, end='')
                print("\n==============================")
                
        except Exception as e:
            print(f"Error during calculation: {str(e)}")
            # If there was an error but the output file exists, show its contents
            if os.path.exists(output_file):
                print("\n=== Output file contents ===")
                with open(output_file, 'r') as f:
                    print(f.read())
                print("===========================")
            raise
        
    except Exception as e:
        print(f"Error during calculation: {str(e)}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
        
    return 0

if __name__ == '__main__':
    main()
