from fabric.api import local, task


@task(default=True)
def test():
    local("nosetests test/**/*.py")


@task
def deploy():
    local("git push origin develop")
    local("git checkout master")
    local("git merge master develop")
    local("git push origin master")
    local("git checkout develop")
