import argparse
from rich_argparse import RichHelpFormatter


class PostGWASHelpFormatter(RichHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


# class PostGWASHelpFormatter(RichHelpFormatter, argparse.RawDescriptionHelpFormatter):
#     """
#     Combined formatter:
#     1. RichHelpFormatter: Renders [colors] and tables.
#     2. RawDescriptionHelpFormatter: Preserves manual newlines (\n).
#     """
#     def __init__(self, *args, **kwargs):
#         # Force a wide width so argparse doesn't wrap lines for you
#         kwargs.setdefault('width', 120)
#         super().__init__(*args, **kwargs)

#     def _fill_text(self, text, width, indent):
#         # This is the 'Nuclear Option': It tells argparse 
#         # "Do NOT touch my text alignment, just return it as is."
#         return "".join(indent + line for line in text.splitlines(keepends=True))



# import argparse
# from rich_argparse import RichHelpFormatter

# class PostGWASHelpFormatter(RichHelpFormatter, argparse.RawDescriptionHelpFormatter):
#     """
#     Forces preservation of newlines and renders Rich colors.
#     """
#     def _fill_text(self, text, width, indent):
#         # Override: Tells argparse to STOP re-wrapping the text.
#         # Just return the lines exactly as they are provided.
#         return "".join(indent + line for line in text.splitlines(keepends=True))