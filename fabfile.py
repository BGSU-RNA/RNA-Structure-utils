from fabric.api import local, task

@task
def test():
    local("nosetests test/**/*.py")
