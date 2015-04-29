import json
import sys

from human_filtering_tools import Blat
from human_filtering_tools import Bwa
from human_filtering_tools import Bmfilter
from human_filtering_tools import Bmtagger
from human_filtering_tools import All_human
from human_filtering_tools import None_human
from human_filtering_tools import Random_human
from human_filtering_tools import Bowtie

""" creates list of tools that can process the input file.
    each tool implements get_human_annotation() method and has name attribute
    and __init__ with one argument.
"""

tools_available = {
    Blat,
    Bwa,
    Bmfilter,
    Bmtagger,
    All_human,
    None_human, 
    Random_human,
    Bowtie
}

def check_if_valid_tool(tool_name, toolname_to_runner):
    if not tool_name in toolname_to_runner:
        raise ValueError("Tool: " + tool_name + 
            " is not a valid tool; change tool name or add it to map of tools(toolname_to_runner)")

def get_parameters_for_tools(path_to_parameter_file="parameters.json"):
    """parse json file with parameters for tools.
    Agrs:
       path_to_params path to json file
    Raises:
        IOError if cannot find or parse file. 
    """
    try:
        params = json.load(open(path_to_parameter_file))
    except ValueError, message:
        print "cannot parse file " + path_to_parameter_file
        print "check for format of the file."
        print "Error message: " + str(message)
        sys.exit(1)
    return params


def create_tools(tool_names, tool_parameters):
    """ create tools for running human filtering.

        For each tool from from the tool_names list checks that tool is valid and create tool
        using tool_parameters dictionary.

        Note: tool_parameters ..
        Returns:
            list of tools(each tool implement runner and parser)

        tool is valid if it is in tools map
        Raises:
             ValueError if tool is not known
    """
    tools = []
    toolname_to_runner = {tool.name : tool for tool in tools_available}
    if not tool_names:
        raise ValueError("empty tools list; need at least one tool")
    for tool_name in tool_names:
        check_if_valid_tool(tool_name, toolname_to_runner)
        try:
            tools.append(toolname_to_runner[tool_name](tool_parameters.get(tool_name, None)))
        except KeyError, message:
            print "Cannot add: " + tool_name + " tool."
            print "Some required parameters are not specified. Check parameters.json file."
            print "Error message: " + str(message)
    return tools
