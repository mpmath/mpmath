import io
import sys


class NullIO(io.BytesIO):

    def __init__(self, *args, **kwargs):
        super(NullIO, self).__init__(*args, **kwargs)

    def write(self, *args, **kwargs):
        pass

    def writelines(self, *args, **kwargs):
        pass


class nostdout(object):

    def __enter__(self):
        self.save_stdout = sys.stdout
        sys.stdout = NullIO()

    def __exit__(self, type, value, tb):
        sys.stdout = self.save_stdout


class stdout(object):

    def __enter__(self):
        pass

    def __exit__(self, type, value, tb):
        pass
