import asyncio
import functools
import subprocess
import os
import numpy as np

class BaseProcessOutputReader:
  """ 
  Class, that run a process, optionally saves all the output of the process to the file
  and pass the stdout and stderr of the process to its two asyn routines, read_error
  and read_output
  """
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
              fd(data)
          fd = stream_reader.feed_data
          stream_reader.feed_data = feed_data

      if self.outfile:
        replace_feed_data(proc.stdout)
        replace_feed_data(proc.stderr)

      out = await asyncio.gather(
          self.read_output(proc.stdout),
          self.read_error(proc.stderr),
          proc.wait(),
        )

      if exception:
         raise exception

      return self.result(*out)

  def result(self, output, error, wait):
      if wait != 0:
         raise ValueError("Process ended with return value {wait}")
      if error is not None:
         return output, error
      return output

  def run(self, cmd, outfile, print_output=False, directory=None, **kwargs):
      self.cmd = cmd
      self.kwargs = kwargs
      self.outfile = outfile
      self.print_output = print_output
      self.directory = directory

      loop = asyncio.get_event_loop()

      import logging
      logging.getLogger("asyncio").setLevel(logging.WARNING)
      out = loop.run_until_complete( self.run_subprocess() )
      if self.outfile:
         self.outfile.close()
      return out

  async def read_error(self, stderr):
      while True:
          line=await stderr.readline()
          if not line:
             return
          print(line.decode('utf8'))

  async def read_output(self, stdout):
      raise NotImplemented('Please, redefine BaseProcess.read_output coroutine')

  def read_from_file(self, output, error=None, return_code=0):

      loop = asyncio.get_event_loop()

      def out():
          air = AsyncioFileReader(output)
          task = loop.create_task(self.read_output(air))
          return loop.run_until_complete(task)

      def err():
          if not error: return None
          air = AsyncioFileReader(error)
          task = loop.create_task(self.read_error(air))
          return loop.run_until_complete(task)

      return self.result(out(), err(), return_code)


class AsyncioFileReader:
  """ `Async' file reader that mimics asyncio StreamReader
      It is in fact not async, but offers the same interface.
  """

  def __init__(self, filename, buffersize=8192):
      self.file = open(filename, 'rb') if isinstance(filename, str) else filename
      self.buffersize = buffersize

  async def readline(self):
      return self.file.readline()

  @staticmethod
  @functools.lru_cache(maxsize=128)
  def separator_shifts(sep):
      shifts = np.ones(len(sep) + 1, dtype=int)
      shift = 1
      for pos in range(len(sep)):
        while shift <= pos and sep[pos] != sep[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift
      return shifts

  async def readuntil(self, sep=b'\n'):
      lsep = len(sep)
      if lsep == 0:
          return b''

      #https://stackoverflow.com/questions/14128763/how-to-find-the-overlap-between-2-sequences-and-return-it
      shifts = self.separator_shifts(sep)

      # do the actual search
      # read the bytes into the array to avoid incrementing the array one by one
      out =  b''
      buffer = bytearray(self.buffersize + lsep)
      bufEnd = 0
      bufPos = 0

      def result():
          if not out:
             return bytes(buffer[:bufPos])
          return out + bytes(buffer[:bufPos])

      startPos = 0
      matchLen = 0
      while True:
        if bufPos == bufEnd:
           if bufEnd >= self.buffersize:
              out += buffer[:bufEnd]
              bufEnd = bufPos = 0
           n = lsep - matchLen
           data = self.file.read(n)
           ldata = len(data)
           buffer[bufPos:bufPos+ldata] = data
           bufEnd += ldata
           if len(data) < n:
              raise asyncio.IncompleteReadError(result(), sep)
        c = buffer[bufPos]
        bufPos += 1
        while matchLen >= 0 and sep[matchLen] != c:
          startPos += shifts[matchLen]
          matchLen -= shifts[matchLen]
        matchLen += 1
        if matchLen == lsep:
           return result()