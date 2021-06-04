from enum import Enum


class FileState(Enum):
    """ Enum that defines the current state of a file."""

    Unread = 1
    Reading = 2
    Read = 3
