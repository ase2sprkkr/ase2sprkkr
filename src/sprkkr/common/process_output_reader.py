import asyncio
import subprocess
import os

class BaseProcessOutputReader:
  """ 
  Class, that run a process, optionally saves all the output of the process to the file
  and pass the stdout and stderr of the process to its two asyn routines, read_error
  and read_output
  """

  def __init__(self, cmd, outfile, print_output=False, directory=None, **kwargs):
      self.cmd = cmd
      self.kwargs = kwargs
      self.outfile = outfile
      self.print_output = print_output
      self.directory = directory

  async def run_subprocess(self):
      loop = asyncio.get_event_loop()

      if self.directory:
         dr = os.getcwd()
         os.chdir(self.directory)
      else:
         dr = None
      proc = await asyncio.create_subprocess_exec(*self.cmd,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        **self.kwargs)
      if dr:
         os.chdir(dr)

      exception = None
      def replace_feed_data(stream_reader):
          nonlocal exception
          if stream_reader._buffer:
             if self.print_output:
                print(stream_reader._buffer.decode('utf8'))
             self.outfile.write(stream_reader._buffer)
          def feed_data(data):
              nonlocal exception
              try:
                if self.print_output:
                  print(data.decode('utf8'))
                self.outfile.write(data)
              except Exception as e:
                if not exception:
                    exception = e
                raise

          stream_reader.feed_data = feed_data

      if self.outfile:
        replace_feed_data(proc.stdout)
        replace_feed_data(proc.stderr)

      await asyncio.gather(
          self.read_output(proc.stdout),
          self.read_error(proc.stderr),
          proc.wait(),
        )
      if exception:
         raise exception



  def run(self):
      loop = asyncio.get_event_loop()

      import logging
      logging.getLogger("asyncio").setLevel(logging.WARNING)
      loop.run_until_complete( self.run_subprocess() )
      if self.outfile:
         self.outfile.close()

  async def read_error(self, stderr):
      while True:
          line=await stderr.readline()
          if not line:
             return
          print(stderr.decode('utf8'))

  async def read_output(self, stdout):
      raise NotImplemented('Please, redefine BaseProcess.read_output coroutine')
