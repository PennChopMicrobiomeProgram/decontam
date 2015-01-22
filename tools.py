from human_filtering_tools import Bmtagger
from human_filtering_tools import All_human
from human_filtering_tools import None_human

""" creates list of tools that can process the input file.
    each tool implements extract_human_reads() method
"""

toolname_to_runner = {
    "bmtagger": Bmtagger(),
    "all_human": All_human(),
    "none_human": None_human()
}


def check_if_valid_tool(tool_name):
    if not tool_name in toolname_to_runner:
        raise ValueError("Tool: " + tool_name + 
            " is not a valid tool; change tool name or add it to map of tools(toolname_to_runner)")


def create_tools(tool_names):
    """
        for each tool from from the tool_names list checks that tool is valid and create tool

        return list of tools(each tool implement runner and parser)

        tool is valid if it is in tools map

        raises ValueError if tool is not known
    """
    tools = []
    if not tool_names:
        raise ValueError("empty tools list; need at least one tool")
    for tool_name in tool_names:
        check_if_valid_tool(tool_name)
        tools.append(toolname_to_runner[tool_name])
    return tools
