import json
import sys

from human_filtering_tools import (
    Snap, Blat, Bwa, Bmfilter, Bmtagger,
    All_human, None_human, Random_human,
    Bowtie,
    )


tools_available = {
    "snap": Snap,
    "blat": Blat,
    "bwa": Bwa,
    "bmfilter": Bmfilter,
    "bmtagger": Bmtagger,
    "all_human": All_human,
    "none_human": None_human,
    "random_human": Random_human,
    "bowtie": Bowtie,
}


def load_parameters(parameter_fp="parameters.json"):
    """parse json file with parameters for tools.
    Agrs:
       path_to_params path to json file
    Raises:
        IOError if cannot find or parse file. 
    """
    return json.load(open(parameter_fp))


def create_tools(tool_names, tool_parameters):
    """Create tool instances from a list of tool names."""
    tools = []
    for tool_name in tool_names:
        tool_cls = tools_available[tool_name]
        params = tool_parameters.get(tool_name, None)
        tool = tool_cls(params)
        tools.append(tool)
    return tools
