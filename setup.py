import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='L-HAPI',
     version='0.1',
     scripts=['L-HAPI.py'] ,
     author="Prajwal Niraula",
     author_email="prajwalniraula@gmail.com",

     description="reduced but fast HAPI version for generating the cross-section",
     long_description=long_description,

     long_description_content_type="text/markdown",
     url="https://github.com/prajwal309/L-HAPI.git",
     packages=setuptools.find_packages(),

     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: MIT License",
         "Operating System :: Tested in LINUX",

     ],

 )
