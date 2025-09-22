"""
This module handles archive creating and uploading to NOMAD
"""

import zipfile
import tempfile
import os
import filecmp
import yaml
from ...outputs.task_result import TaskResult
from ...common.decorators import cached_property
from ...common.yaml import IndentDumper
from typing import Optional, Union, Dict


def map_io_to_nomad(items):

    return [ {'name': i[0], 'section' : i[1] }
             if isinstance(i, tuple) else {'section' : i }
             for i in items ]


class ExternalEntry():

    def __init__(self, outputs):

        def make_full(v):
            if '#' in v:
                return v
            return v + "#/data/outputs[0]"

        if not isinstance(outputs, list):
            self.outputs = [ make_full(outputs) ]
        else:
            self.outputs = [ make_full(o) for o in outputs ]

    def outputs(self):
        yield from self.outputs


class NomadEntry():
    """ Description of a entry in a :class:`NomadArchive` """

    def __init__(self, archive, name, output:TaskResult, depends_on:"Union[str, bool, NomadEntry]"):
        self.archive = archive
        self.name = name
        self.output = output
        self.depends_on = depends_on

    @property
    def depends_on(self):
        return self._depends_on

    @depends_on.setter
    def depends_on(self, value):
        if isinstance(value, NomadEntry):
            p1 = self.output.path_to('potential')
            p2 = value.output.path_to('converged')
            if p1 != p2 and not filecmp.cmp(p1, p2):
               raise ValueError(f"The task {self.output.files['output']} does not use "
                                f"potential from {value.output.files['output']}")
            self.archive.add_symlink(os.path.join('..', value.file('converged')), self.file('potential'))
        self._depends_on = value

    def file(self, name):
        return os.path.join(os.path.dirname(self.name), os.path.basename(self.output.files[name]))

    @cached_property
    def task_name(self):
        return self.output.task_name.upper()

    @cached_property
    def symbols(self):
        return str(self.output.potential.atoms.symbols)

    def _resource(self, resource):
        return f'../upload/archive/mainfile/{self.name}#/{resource}/'

    def outputs(self):
        yield f"{self.task_name} of {self.symbols}", self._resource('data/outputs[0]')

    def inputs(self):
        yield self.symbols, self._resource('data/model_system[0]')
        yield self.task_name, self._resource('data/model_method[0]')
        if self.depends_on:
            yield from self.depends_on.outputs()

    def task(self):
        return {
            'name'  : f"{self.output.task_name} for {self.symbols}",
            'm_def' : 'nomad.datamodel.metainfo.workflow.TaskReference',
            'inputs': map_io_to_nomad(self.inputs()),
            'outputs':map_io_to_nomad(self.outputs()),
        }


class NomadArchive():
    """ This class handles Nomad uploads """

    def __init__(self, filename:Optional[str]=None, depends:Union[str, bool]=True, name=None):
        """
        Parameters
        ----------
        filename
          Name of the resulting zip archive

        depends
          Added tasks will be (by default) dependendent on a given entry.
          False means no dependency
          True means autodetect -- only one SCF task can be added and this will be the dependency

        name
          Name of the whole workflow
        """
        if filename is None:
            self.file = tempfile.TemporaryFile(mode='w+b')
        else:
            self.file = open(filename, 'w+b')
        self.zip = zipfile.ZipFile(self.file, 'x', zipfile.ZIP_BZIP2)
        self.entries = {}
        self.scf = None
        self.depends = ExternalEntry(depends) if isinstance(depends, str) else depends
        self.name = name

    def add_symlink(self, frm, to):
        zi = zipfile.ZipInfo(to)
        zi.create_system = 3
        zi.external_attr = 2716663808
        self.zip.writestr(zi, frm)

    def _add_entry(self, output, depends):

        if not output.files['output']:
            raise ValueError('Output file name has to be specified')

        file = os.path.basename(output.files['output'])
        folder = folder_base = os.path.splitext(file)[0]
        counter = 1
        while folder in self.entries:
            folder = f"{folder_base}_{self.counter}"
            counter+=1
        self.zip.mkdir(folder)
        name = f"{folder}/{file}"

        if isinstance(depends, str):
            depends = ExternalEntry(depends)

        out = NomadEntry(self, name, output, depends)

        def add_file(fname):
            file = os.path.basename(fname)
            self.zip.write(output.path_to(kind), os.path.join(
              folder, file
            ))

        for kind in output.files:
            if kind == 'potential' and \
               not isinstance(depends, ExternalEntry) and \
               not output.task_name.lower() == 'scf':
                  continue
            add_file(output.path_to(kind))

        self.entries[name] = out
        return out

    def add_entry(self, output:Union[TaskResult,str], depends:Union[str,bool, NomadEntry]=True):
        """
        Add entry

        Parameters
        ----------
        output
          Output file to add

        depends
          ``str``:  "foreign entry point"
          NomadEntry: Already added package
          ``True``:  Automatic detection
          ``False``: No dependency
        """
        if not isinstance(output, TaskResult):
            output = TaskResult.from_file(output)

        if depends is True:
            if output.task_name.lower() == 'scf':
                depends = False
            else:
                depends = self.depends

        if depends is False and output.task_name.lower() != 'scf':
            raise ValueError("Non-SCF-task have to depend on some SCF task")
        out = self._add_entry(output, depends)
        return out

    def finalize(self):
        self.resolve_auto_dependencies()
        self.zip.writestr('workflow.archive.yaml', yaml.dump(self.workflow(),
                                                   Dumper=IndentDumper,
                                                   sort_keys=False))
        self.zip.close()
        self.file.seek(0)

    def resolve_auto_dependencies(self):
        """ If there is any entry with 'auto' dependency,
        make it dependent to a SCF task, which have to be uniqe
        """
        auto = [ i for i in self.entries.values() if i.depends_on is True ]
        if auto:
            scf = [ i for i in self.entries.values() if i.output.task_name.lower() == 'scf' ]
            if len(scf) != 1:
                raise ValueError("If an automatic dependency is used, just one SCF task have to be submited.")
            scf = scf[0]
            for i in auto:
                i.depends_on = scf

    def workflow(self)->Dict:
        """
        Returns
        -------

        workflow
          Dictionary describing YAML for NOMAD workflow
        """
        def gather(what):
            out = set()
            for i in self.entries.values():
                out |= set(getattr(i, what)())
            return list(out)

        def tasks():
            return [
              i.task() for i in self.entries.values()
            ]

        workflow = {
            'm_def': 'nomad.datamodel.metainfo.workflow.TaskReference',
            'inputs':map_io_to_nomad(gather('inputs')),
            'outputs':map_io_to_nomad(gather('outputs')),
            'tasks': tasks()
        }

        out = { 'workflow2' : workflow }
        if self.name:
            out['name'] = self.name
        return out
