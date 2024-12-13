---
sort: 2
---

This list has been generated using `git shortlog -s --all`.
If your name is missing, please let philip.cardiff@ucd.ie know.

* Philip Cardiff
* Zeljko Tukovic
* Danial Khazaei
* Denis Maier
* Iago Lessa de Oliveira
* Jozsef Nagy
* Jeff Heylmun
* Saber Mohammadi
* Bruno Santos
* Cyrille Bonamy
* Jeffrey Heylmun
* Scott Levie
* Sindeev Sergey
* Andrew Whelan
* Emad Tandis
* Karen FitzGerald
* Philipp Milovic
* Seevani Bali
* Ivan Batistic
* Simona Moretti
* Andrea Luigi Facci
* Federico Mazzanti
* Xiaohu Guo
* Amirhossein Taran
* Mark Olesen
* Sandro Manservisi
* Greg Burgreen

# Contributing to solids4foam

Thank you for considering contributing to **solids4foam**!
Thanks for helping create a world where fluid and structure are solved using your favorite CFD code.

## How Can I Contribute?

### 1. Reporting Bugs

If you encounter a bug, please open an [GitHub issue](https://github.com/solids4foam/solids4foam/issues) and include:

- A concise and descriptive title.

- Steps to reproduce the issue along with and explanation of the expected behavior.

- Relevant details, such as screenshots, logs, or error messages and an environment description in the following form:

    ```plaintext
    - OS version: Ubuntu 22.04.5 LTS
    - OpenFOAM version: openfoam2312
    - solids4foam version: 2.1
    - Docker image (if applicable)
    ```

### 2. Suggesting Features

We welcome ideas for new features! To suggest a feature:

- Open a [GitHub issue](https://github.com/solids4foam/solids4foam/issues) with title: **Feature request: Feature title**.
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

#### Step 2: Create a Branch

Create a new branch for your changes:

```bash
git checkout -b branch_name
```

Branch naming conventions:

- **New features**: Start with `feature-name` (e.g., `feature-newSolver`).
- **Bug fixes**: Start with `fix-name` (e.g., `fix-boundaryConditions`).
- **Code refactoring**: Start with `refactor-name` (e.g., `refactor-meshHandling`).
- **Documentation changes**: Start with `docs-name` (e.g., `docs-update-readme` or `docs-fix-typos`).

#### Step 3: Make Your Changes

- Follow the [OpenFOAM Coding Style Guide](https://openfoam.org/dev/coding-style-guide/)).
- Ensure that your code can be compiled and adheres the existing code structure and organization.
- Test your changes thoroughly with relevant test cases. Include or update tutorial cases if applicable.
- Document your code with clear comments.
- Before submitting, check the code using available linters or style checkers ([shellcheck](https://www.shellcheck.net/#), [markdownlint](https://dlaa.me/markdownlint/)).

#### Step 4: Commit Your Changes

Write clear and descriptive commit messages and commit your changes. Reference any related issues using hashtags if applicable::

```bash
git add file_name.C
git commit -m "Fix for issue #51"
```

#### Step 5: Push and Open a Pull Request

Push your branch to your fork:

```bash
git push origin branch_name
```

Go to the original solids4foam repository on GitHub and open a pull request (target the development branch for your pull request). Include:

- **A clear title**, such as `Fix: Issue #123` or `Feature: Add new functionality`.
- **A detailed description** of your changes.
- Any **related issue numbers** (e.g., `Closes #123`).
- **Instructions or example cases** to test your contribution, if applicable.

## Questions?

If you have any questions or need assistance, feel free to open a [discussion](https://github.com/solids4foam/solids4foam/discussions).

Thank you for contributing to solids4foam!
