import sys

def main():
    from ase2sprkkr.sprkkr.calculator import SPRKKR
    import os

    print("Starting SPRKKR calculation...")

    # Initialize calculator
    calculator = SPRKKR()

    # Set input parameters
    print("Setting up input parameters...")
    calculator.input_parameters = 'xas'
    calculator.input_parameters.CONTROL.DATASET = 'Fe'
    calculator.input_parameters.TASK.IT = 1  
    calculator.input_parameters.TASK.FRAMETET= 90.0
    calculator.input_parameters.TASK.FRAMEPHI = 45.0
    current_dir = os.getcwd()

    try:
        potential_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Fe.pot_new')
        result = calculator.calculate(
            potential=potential_file,
            print_output=True,
            output_file=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Fe_XAS.out')
        )

    except Exception as e:
        print(f"Error during calculation: {str(e)}")
        raise

    return 0

if __name__ == '__main__':
    main()
