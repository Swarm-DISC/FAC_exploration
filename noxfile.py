"""
For example:
nox --session lint
nox --session build
"""
import nox


@nox.session
def lint(session: nox.Session) -> None:
    """
    Run the linter.
    """
    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files", *session.posargs)


@nox.session(venv_backend="mamba")
def build(session: nox.Session) -> None:
    """
    Generate the Jupyter Book
    """
    session.run(
        *"mamba env update --file environment.yml".split(" ")
    )
    session.run("bash", "build_book.sh")


