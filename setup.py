from setuptools import setup, find_packages
from chime import version

with open('README.md', 'r') as fh:
    long_description = fh.read()

name = 'chime'

setup(
    name=name,
    version=version.CHIME_VERSION,
    author='Ryuichiro Hara',
    author_email='hara.ryuichiro@gmail.com',
    description=name,
    long_description=long_description,
    long_description_content_type='text/markdown',
    packages=find_packages(),
    scripts=['tools/atomic_conv.py', 'tools/chime_test.py'],
    url='https://gitlab.com/rhara/%s' % name,
    classifiers = [
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)
