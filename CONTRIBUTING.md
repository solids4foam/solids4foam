# Contributing to solids4foam

Thank you for considering contributing to **solids4foam**!

Thanks for helping create a world where fluid and structure are solved using your favourite CFD code.

## How Can I Contribute?

### 1. Reporting Bugs

If you encounter a bug, please open a [GitHub issue](https://github.com/solids4foam/solids4foam/issues) and include the following:

- A concise and descriptive title.

- Steps to reproduce the issue, along with explaining the expected behaviour.

- Relevant details, such as screenshots, logs, or error messages and an environment description in the following form:

    ```plaintext
    - OS version: Ubuntu 22.04.5 LTS
    - OpenFOAM version: openfoam2312
    - solids4foam version: 2.1
    - Docker image (if applicable)
    ```

### 2. Suggesting Features

We welcome ideas for new features! To suggest a feature:

- Open a [GitHub issue](https://github.com/solids4foam/solids4foam/issues) with the title: **Feature request: feature title**.
- Provide a detailed description of the proposed feature.
- Explain its potential benefits for the solids4foam toolbox.

### 3. Submitting Code Changes

#### Step 1: Fork and Clone

1. Fork the repository by clicking the "Fork" button on GitHub.

2. Clone your fork locally:

    ```bash
    git clone https://github.com/your-username/solids4foam.git
    cd solids4foam
    ```

#### Step 2: Make Your Changes

- Follow the [OpenFOAM Coding Style Guide](https://openfoam.org/dev/coding-style-guide/)).
- Ensure your code can be compiled and adheres to the existing code structure and organisation.
- Test your changes thoroughly with relevant test cases. Include or update tutorial cases if applicable.
- Document your code with clear comments.
- Before submitting, check the code using available linters or style checkers ([shellcheck](https://www.shellcheck.net/#), [markdownlint](https://dlaa.me/markdownlint/)).

#### Step 3: Commit Your Changes

Write clear and descriptive commit messages and commit your changes. Reference any related issues using hashtags if applicable::

```bash
git add file_name.C
git commit -m "Fix for issue #51"
```

#### Step 4: Push and Open a Pull Request

Push your branch to your fork:

```bash
git push
```

Go to the original solids4foam repository on GitHub and open a pull request (target the `development` branch for your pull request). Include:

- **A clear title**, such as `Fix: Issue #123` or `Feature: Add new functionality`.
- **A detailed description** of your changes.
- Any **related issue numbers** (e.g., `Closes #123`).
- **Instructions or example cases** to test your contribution, if applicable.

## Questions?

If you have any questions or need assistance, we encourage you to use the appropriate discussion platform based on your query:

- **User Support**: For general usage questions, user support, or troubleshooting, please post on the **CFD Online Forum**: [https://www.cfd-online.com/Forums/openfoam-cc-toolkits-fluid-structure-interaction/](https://www.cfd-online.com/Forums/openfoam-cc-toolkits-fluid-structure-interaction/)

- **Developer Discussions**: For discussions about development, contributions, or feature requests, please open a discussion on our **GitHub Discussions** page: https://github.com/solids4foam/solids4foam/discussions](https://github.com/solids4foam/solids4foam/discussions)

Thank you for contributing to solids4foam!
