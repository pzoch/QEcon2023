# QEcon2023


> [!IMPORTANT]
> **Problem Set 2** is available. It is due on **18.12.2023, 11:59 PM**. You can submit your code by email or (preferably) by sharing a self-contained repository with all required files inside.

## Log and annoucements 
+ 4/12: we proved the equivalence between solving the sequence problem and the functional (Bellman) equation for (finite) Markov Decision Processes. We introduced policy operators and talked about Howard policy iteration and optimistic policy iteration. 
+ 27/11: two examples of dynamic programming: option pricing and inventory management problem.
+ 20/11: we finished intro to dynamic programming by using value function iteration to solve the McCall's job search model and talked about Markov Chains.
+ 13/11: we started discussing (discrete) dynamic programming. We introduced the McCall's job search model with finite and infinite horizon and solved it by iterating on the continuation value. Then we talked about the Contraction Mapping Theorem and how it potentially allows us to solve Bellman equations.
+ 7/11: we talked a bit more about spectral methods (mainly Chebyshev polynomials), we saw some interpolation code example and finished by discussing the finite element method.
+ 23/10: we finished the topic of optimization by spending some time on constrained optimization; we saw how to draw random numbers in Julia and started talking about spectral methods.
+ 16/10: we talked about multidimensional unconstrained optimization: Newton, polytope methods, direction methods.
+ 9/10: we talked about unidimensional unconstrained optimization: bracketing methods and Newton's method.

## How to start using Julia
Follow [these instructions](https://code.visualstudio.com/docs/languages/julia). Please use Julia 1.9.3.

Once you completed all these steps press <kbd>Alt</kbd>+<kbd>j</kbd> <kbd>Alt</kbd>+<kbd>o</kbd> in VS Code (alternatively press <kbd>Ctrl</kbd>+<kbd>Shift</kbd>+<kbd>p</kbd> to get the command palette and find the shortcut). This should start Julia REPL and you will see green `julia` text in the terminal. Input <kbd>]</kbd> (right square bracket) in the terminal. This should change text from green to blue (and you should see something like `1.9`). Type `add IJulia` and press <kbd>enter</kbd>. This will download a package that allows you to use Jupyter notebooks with Julia.

### Environments
[This note](https://jkrumbiegel.com/pages/2022-08-26-pkg-introduction/) might be helpful.

## Git tutorial 
Please see [this tutorial](https://swcarpentry.github.io/git-novice/).   

