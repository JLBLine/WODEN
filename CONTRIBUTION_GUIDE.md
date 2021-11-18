# Contributing to WODEN

If you have a cool feature request, a burning passion to add some functionality
to `WODEN`, or have found an annoying bug, please reach out. I welcome the help!
Here is a quick guide on how to go about it.

## I have an idea for a feature
Please raise an issue suggesting the new feature, and chuck the `enhancement`
label on the issue. Try and describe the feature as completely as possible and
why it would be useful for you. We can then start a dialogue and work out the
best way forward / where in the priority list of new features to implement this
it might lie.

## I found a bug
`#sadface`, these things happen. Please raise an issue and add the `bug` label on the issue (please check the issue doesn't already exist). Include a description of what's going wrong, including an example command that causes the error for you.

## I have an idea/found a bug AND I want to fix it myself
Awesome, go for it! First off, make sure you've submitted an issue so we know
the bug/feature isn't being fixed or implemented in another branch somewhere. After
a little discussion to clarify the goals, make your own branch or fork, and
then start working on that branch/fork. For any piece of work, I'd ask that:

 - If fixing code, ensure the fix still passes any existing unit test as detailed in [the online documentation here](https://woden.readthedocs.io/en/joss_review/testing/cmake_testing.html#what-do-the-tests-actually-do). If the fix changes the expected outcomes, discuss it on the issue, and if appropriate, you can update the existing unit tests.
 - If making new code, write new unit tests that cover as much of the new code as possible.
 - Document whatever you've done so we all know what good stuff you've added in.

If you've not written test code before, or made documentation via `sphinx`, just let me know on whatever git issue you've started (maybe use the "help wanted" label) and we can work together on it.

Once all that's done, submit a pull request and we'll get it in the main branch.

## Response time
At the moment, it's just me (Jack Line) writing the code, so please
be patient with me and I'll get back to you as fast as my schedule allows.
