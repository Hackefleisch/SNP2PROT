"""Model interface. TODO."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class BaseModel(ABC):
    """Minimal contract every model in this project satisfies."""

    @abstractmethod
    def fit(self, train_data: Any, val_data: Any | None = None) -> BaseModel: ...

    @abstractmethod
    def predict(self, data: Any) -> Any: ...
