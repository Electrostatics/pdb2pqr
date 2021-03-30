# Contributing to PDB2PQR

We are excited to have you contribute to PDB2PQR!
If you have general questions about the code, please see the [documentation](https://pdb2pqr.readthedocs.io/) and/or start a thread in the [discussion forum](https://github.com/Electrostatics/pdb2pqr/discussions). 
Once you have familiarity with the code, you're ready to start contributing:

## Find or create a problem

Find a problem/feature that needs to be resolved in the code [issues](https://github.com/Electrostatics/pdb2pqr/issues).
If the problem you're trying to solve doesn't exist, create a Github issue to resolve some or all of the problem or to add features.
For example, you could
  * Fix a bug
  * Add a new capability with 3 parts - create issue for each part it possible
  * In general, make problem/feature small enough to finish it within a week or two

## Create a branch

Create a git branch using the convention of `github_user`/`issue_#` where `github_user` is your user name and `issue_#` is the issue number from step 2 above.

## Create and pass tests

Create a test that replicates the problem/feature and fails and show how your fix results in a working test.
Work on the code until your test, as well as all previous tests, pass.

## Do not submit messy code

Run your code through the following steps:

1. Formatting tools:  [psf/black](https://github.com/psf/black)
2. Linting tools: [pylint](https://www.pylint.org/) and [flake8](https://flake8.pycqa.org/en/latest/)

## Submit a pull request

1. Commit your changes to Git and push your branch to the Github repo.
2. Create a [pull request](https://github.com/Electrostatics/pdb2pqr/compare?expand=1) and add reviewers (e.g., `speleo3`, `orbeckst`, `sobolevnrm`, and `intendo`) to the request.
[Reference the issue you are trying to fix](https://docs.github.com/en/github/managing-your-work-on-github/linking-a-pull-request-to-an-issue) in your pull request.
3. If the pull request passes [Github Actions](https://github.com/features/actions) and peer review, then the branch will be merged with `main`, your branch will be deleted, and the development team will be very grateful for your contribution!

Thank you for considering to contribute to our code!
