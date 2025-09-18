### Contribute to the Software

There are many ways to contribute to ggoncoplot.

1.  Request features you would like to by [creating new issues on
    github](https://github.com/selkamand/ggoncoplot/issues)
2.  [Make your visualisation packages ggoncoplot-compatible]
3.  [Directly contribute to the ggoncoplot codebase]

#### Make your visualisation packages ggoncoplot-compatible

If your package produces ggplots that you would like to interactively
link with ggoncoplot, consider converting your geoms to their ggiraph
interactive equivalents and adding a data_id based on a sample
identifier column in the dataset. That way end-users can create a
data-linked oncoplot composed with your packages plots using patchwork
([example](https://selkamand.github.io/ggoncoplot/articles/manual.html#interaction-with-other-packages))

#### Directly contribute to the ggoncoplot codebase

We welcome contributions from the community to enhance and expand the
functionality of `ggoncoplot`. Whether you want to fix a bug, add new
features, improve documentation, or optimize performance, your efforts
are highly valued. To get started:

1.  **Fork the Repository**: Click on the 'Fork' button at the top right
    of this page to create a copy of the repository in your GitHub
    account.

2.  **Clone the Repository**: Use `git clone` to clone your forked
    repository to your local machine.

    ``` bash
    git clone https://github.com/selkamand/ggoncoplot.git
    ```

3.  **Create a Branch**

    ``` bash
    git checkout -b feature-name
    ```

4.  **Make Changes**: Implement your changes in the new branch

5. **Test Package**: Run `devtools::test()` and `devtools::check()` to ensure changes do not break other behaviours. 
Add additional unit tests to cover new behaviour.

7.  **Commit and Push**: Commit your changes and push the branch to your
    forked repository.

8.  **Create a Pull Request**: Go to the original repository and open a
    pull request from your branch. Please provide a clear description of
    your changes and any relevant issues or discussions.

### Report Issues or Problems with the Software

If you encounter any issues, bugs, or have suggestions for improvements,
please report them using the [GitHub Issues Tab](https://github.com/selkamand/ggoncoplot/issues/).

### Seek Support

For any questions or support regarding the use of ggoncoplot you can:

-   **Check the Documentation**: Comprehensive documentation is
    available [here](https://selkamand.github.io/ggoncoplot/index.html).

-   **Create a** [new
    issue](https://github.com/selkamand/ggoncoplot/issues/new) with your
    query.

-   **Browse Existing Issues**: Check the
    [Issues](https://github.com/selkamand/ggoncoplot/issues) page to see
    if your query has been addressed.

-   **Contact Us**: If you need direct assistance, please [contact the
    maintainers
    directly](mailto:selkamand@ccia.org.au?subject=ggoncoplot)
