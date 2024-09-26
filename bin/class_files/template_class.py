# Author: Moritz Reichert
# Date  : 26.09.2024
import numpy as np


class template(object):
    """
      Class to read a WinNet template file.
    """

    def __init__(self, path):
        """
          Initialize the template class.
        """
        self.path = path


    def read_data(self):
        """
          Read the data from the template file and store it in a dictionary.
        """
        # Create an empty dictionary to store the entries
        self.__entries = {}

        # Read the data from the file
        with open(self.path, 'r') as f:
            self.data = f.readlines()
            for line in self.data:
                if line.strip().startswith('#'):
                    continue

                if line.strip() =="":
                    continue

                key = line.split("=")[0].strip()
                value = line.split("=")[1].strip().replace('"', '')
                self.__entries[key] = value

    @property
    def entries(self):
        """
          Get the entries of the template file.
        """
        if not hasattr(self, '__entries'):
            self.read_data()
        return self.__entries


    def __getitem__(self, key):
        """
          Get the value of a specific key.
        """
        if not hasattr(self, 'entries'):
            self.read_data()
        return self.entries[key]


if __name__ == '__main__':
    # Example:
    path = '../../runs/test'
    t = template(path)
    print(t["isotopes_file"])