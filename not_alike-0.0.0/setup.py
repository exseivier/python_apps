#!/usr/bin/env python3.9


from setuptools import setup

setup(
    name='not_alike',
    version='0.0.0',
    author='Javier Montalvo',
    author_email='buitrejma@gmail.com',
    py_modules=['man_cli'],
    python_requires='>=3.9',
    description='This package is an interprocess manager of a pipeline that finds not alike regions of query genome compared to a hugh list of different genomes',
    install_requires =['click'],
    scripts = ['scripts/not_alike.sh'],
    url='https://www.github.com/exseivier/',
    entry_points = '''
    [console_scripts]
    manager=man_cli:man_cli
    '''
        )
