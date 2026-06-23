from pathlib import Path

from pydantic import BaseModel


class Settings(BaseModel):
    project_root: Path = Path.cwd()
    data_root: Path = Path("data")
    output_root: Path = Path("outputs")
    config_root: Path = Path("configs")

