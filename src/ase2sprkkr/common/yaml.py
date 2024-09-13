""" YAML - related stuff """
import yaml


class IndentDumper(yaml.Dumper):
    """
    Better indentation for YAML, according to

    https://stackoverflow.com/questions/25108581/python-yaml-dump-bad-indentation
    """
    def increase_indent(self, flow=False, indentless=False):
        return super().increase_indent(flow, False)
