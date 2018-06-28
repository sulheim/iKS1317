# Make model file in all formats

import cobra
import os
import logging
from subprocess import check_output
import sys

def get_name_without_extension(path):
    folder_path, name = os.path.split(path)
    name, ext = os.path.splitext(name)
    return name


def get_path(sbml_path, new_extension, new_name = None):
    folder_path, name = os.path.split(sbml_path)
    name, ext = os.path.splitext(name)

    if isinstance(new_name, str):
        name = new_name

    return os.path.join(folder_path, name+new_extension)

def run(sbml_path, logger = None):
    if not logger:
        logger = logging.getLogger(__name__)

    logger.info("Reading SBML '%s'.", sbml_path)
    model = cobra.io.read_sbml_model(sbml_path)
    json_path = write_json(model, sbml_path, logger)
    matlab_path = write_matlab(model, sbml_path, logger)
    legacy_path = write_legacy_cobra(model, sbml_path, logger)
    no_fbc_path = write_sbml_no_fbc(model, sbml_path, logger)

    
    check_output(["git", "add", json_path])
    check_output(["git", "add", matlab_path])
    check_output(["git", "add", legacy_path])
    check_output(["git", "add", no_fbc_path])
    

def write_json(model, sbml_path, logger):
    json_path = get_path(sbml_path, ".json")
    logger.info("Writing json '%s'.", json_path)
    cobra.io.save_json_model(model, json_path)
    return json_path

def write_matlab(model, sbml_path, logger):
    matlab_path = get_path(sbml_path, ".mat")
    logger.info("Writing matlab '%s'.", matlab_path)
    cobra.io.save_matlab_model(model, matlab_path)
    return matlab_path

def write_legacy_cobra(model, sbml_path, logger):
    new_name = get_name_without_extension(sbml_path) +  "_legacy"
    new_path = get_path(sbml_path, ".xml", new_name)
    logger.info("Writing legacy cobra '%s'.", new_path)
    cobra.io.write_legacy_sbml(model, new_path)
    return new_path
    
def write_sbml_no_fbc(model, sbml_path, logger):
    new_name = get_name_without_extension(sbml_path) +  "_no_fbc"
    new_path = get_path(sbml_path, ".xml", new_name)
    logger.info("Writing cobra without fbc '%s'.", new_path)
    cobra.io.write_sbml_model(model, new_path, use_fbc_package = False)
    return new_path





if __name__ == '__main__':
    if len(sys.argv) == 1:
        path = r"C:/Users/snorres/git/gem_sco/iKS1317.xml"
    else:
        path = sys.argv[1]

    run(path)