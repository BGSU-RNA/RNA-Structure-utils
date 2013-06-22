from fabric.api import local, task


@task(default=True)
def test():
    local("nosetests")


@task
def tox():
    local("tox")


@task
def deploy():
    tox()
    local("git push origin develop")
    local("git checkout master")
    local("git merge master develop")
    local("git push origin master")
    local("git checkout develop")
