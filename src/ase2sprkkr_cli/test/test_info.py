from ase2sprkkr_cli.commands import info as info_command
import argparse

def test_info_lists_stable_tasks_before_experimental(capsys):
    parser = argparse.ArgumentParser()
    info_command.parser(parser)
    info_command.run(parser.parse_args(["--task"]), None)
    output = capsys.readouterr().out
    heading = "  Experimental tasks\n  ~~~~~~~~~~~~~~~~~~"
    assert output.count(heading) == 1
    regular, experimental = output.split(heading)
    assert "SCF" in regular
    assert "BSFEK" in regular
    assert "BSFKK" in regular
    assert "AES" in experimental
    assert "DMI" in experimental
    assert "PLOTSFN" in experimental
    assert " (experimental)" not in output

    info_command.run(parser.parse_args(["--task=dmi"]), None)
    filtered = capsys.readouterr().out
    assert heading in filtered
    assert "DMI" in filtered
