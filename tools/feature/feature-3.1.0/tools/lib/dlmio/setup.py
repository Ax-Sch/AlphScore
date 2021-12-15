from distutils.core import setup,Extension

module1 = Extension('dlmio',sources=['dlmio.c'])

setup(name = 'dlmio',
        version = '1.0', 
        description = 'Delimited file reading', 
        ext_modules = [module1])
