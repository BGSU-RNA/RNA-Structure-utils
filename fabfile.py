from fabric.api import local, task

@task(default=True)
def test():
    local("nosetests test/**/*.py")
