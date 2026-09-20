"""Post-hoc analysis of pipeline output.

Reads the tables the pipeline writes; nothing here runs as a pipeline step. Its
plotting and stats dependencies live in the ``analysis`` environment, so
``pixi install -e analysis`` is needed before these modules import.
"""
