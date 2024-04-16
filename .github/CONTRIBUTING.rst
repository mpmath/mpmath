Reporting Issues
================

When opening a new issue, please take the following steps:

1. Please search `GitHub issues`_ to avoid duplicate reports.

2. If possible, try updating to master and reproducing your issue.

3. Try to include a minimal code example that demonstrates the problem.

4. Include any relevant details of your local setup (mpmath version, Python
   version, installed libraries).

Please avoid changing your messages on the GitHub, unless you want fix a typo
and so on.  Just expand your comment or add a new one.


Contributing Code
=================

All work should be submitted via `Pull Requests (PR)`_.

1. PR can be submitted as soon as there is code worth discussing.
   Please make a draft PR, if one is not intended to be merged
   in its present shape even if all checks pass.

2. Please put your work on the branch of your fork, not in the
   master branch.  PR should generally be made against master.

3. One logical change per commit.  Make good commit messages: short
   (<= 78 characters) one-line summary, then newline followed by
   verbose description of your changes.  Please `mention closed
   issues`_ with commit message.

4. Please conform to `PEP 8`_; run::

       flake518

   to check formatting.

5. PR should include tests:

   1. Bugfixes should include regression tests (named as ``test_issue_123``).
   2. All new functionality should be tested, every new line
      should be covered by tests.  Please use in tests only
      public interfaces.  Regression tests are not accounted in
      the coverage statistics.
   3. Optionally, provide doctests to illustrate usage.  But keep in
      mind, doctests are not tests.  Think of them as examples that
      happen to be tested.

6. It's good idea to be sure that **all** existing tests
   pass and you don't break anything, so please run::

       pytest

7. If your change affects documentation, please build it by::

       sphinx-build -W -b html docs build/sphinx/html

   and check that it looks as expected.


.. _GitHub issues: https://github.com/mpmath/mpmath/issues
.. _Pull Requests (PR): https://github.com/mpmath/mpmath/pulls
.. _PEP 8: https://www.python.org/dev/peps/pep-0008/
.. _mention closed issues: https://help.github.com/en/github/managing-your-work-on-github/linking-a-pull-request-to-an-issue
