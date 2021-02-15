from sprkkr.calcio import (InputFile,
                           load_parameters, make_sections, compare_sections)

test_file = InputFile('rw-test.inp', defaults_filename='defaults/scf.inp')
test_file.control_section.set(DATASET="Fe", ADSI="SCF", POTFIL="Fe.pot")
test_file.write()

pars = load_parameters('rw-test.inp')
sections = make_sections(pars)

print(sections['control'] == test_file.sections['control'])

print(compare_sections(sections, test_file.sections))
