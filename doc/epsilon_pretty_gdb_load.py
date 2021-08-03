import gdb.printing
import epsilon_pretty_printing

gdb.printing.register_pretty_printer(
    gdb.current_objfile(),
    epsilon_pretty_printing.build_pretty_printer())