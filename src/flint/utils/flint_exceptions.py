

class DomainError(Exception):
    """
    Exception intended to be called when a method is called on a
    ring or field for which the domain is invalid.
    """
    pass


class IncompatibleContextError(Exception):
    """
    Exception intended to be called when a method involves two or more
    incompatible contexts.
    """
    pass
