
class CompilationError(Exception):
    """Represents a single error found during compilation (lexing, parsing, or semantic analysis)."""
    def __init__(self, message, line=None, column=None):
        super().__init__(message)
        self.message = message
        self.line = line
        self.column = column

    def __str__(self):
        if self.line and self.column:
            return f"[Line {self.line}, Column {self.column}] {self.message}"
        elif self.line:
            return f"[Line {self.line}] {self.message}"
        return self.message

class ErrorReporter:
    """A class to collect and report compilation errors."""
    def __init__(self):
        self._errors = []

    def has_errors(self):
        """Returns True if any errors have been reported."""
        return len(self._errors) > 0

    def report(self, message, line=None, column=None):
        """Records a new compilation error."""
        self._errors.append(CompilationError(message, line, column))

    def get_errors(self):
        """Returns the list of collected errors."""
        return self._errors

    def display(self):
        """Prints all collected errors to stderr."""
        import sys
        if not self.has_errors():
            return
        
        print(f"Compilation failed with {len(self._errors)} error(s):", file=sys.stderr)
        for error in sorted(self._errors, key=lambda e: (e.line or 0, e.column or 0)):
            print(f"- {error}", file=sys.stderr)

    def __iter__(self):
        return iter(self._errors)
