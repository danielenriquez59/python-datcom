"""Configuration management for PyDATCOM."""

from pathlib import Path
import yaml

_config_dir = Path(__file__).parent


def load_defaults():
    """Load default configuration from YAML file."""
    defaults_path = _config_dir / 'defaults.yaml'
    with open(defaults_path, 'r') as f:
        return yaml.safe_load(f)


__all__ = ['load_defaults']

