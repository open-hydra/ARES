# Changelog

All notable changes to ARES are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- A [Versioning](docs/development/versioning.md) page: the releases ARES requires of its libraries, and how a build gets them.

### Changed

- FiNeR moved to `lib/FiNeR`, pinned at v2.3.0; its libraries are fetched by the configure, no longer submodules.
- A library at another release than ARES requires stops the configure (`-DARES_PIN_CHECK=WARNING` to override).

## [0.1.0-alpha]

First tagged version.
