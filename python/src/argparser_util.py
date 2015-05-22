"""
Created on May 22, 2015

@author: Steven Wu
"""
import argparse


class CustomSingleMetavarFormatter(argparse.HelpFormatter):

    def __init__(self, prog,
                 indent_increment=2,
                 max_help_position=50,
                 width=None):
        super(CustomSingleMetavarFormatter, self).__init__(prog, indent_increment, max_help_position, width)

#     self._max_help_position = 50
    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            # if the Optional doesn't take a value, format is:
            #    -s, --long
            if action.nargs == 0:
                parts.extend(action.option_strings)

            # if the Optional takes a value, format is:
            #    -s ARGS, --long ARGS
            # change to
            #    -s, --long ARGS
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    # parts.append('%s %s' % (option_string, args_string))
                    parts.append('%s' % option_string)
                parts[-1] += ' %s' % args_string
            return ', '.join(parts)
