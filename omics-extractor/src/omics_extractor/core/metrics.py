from __future__ import annotations
from typing import Any, Dict, Literal, Optional
from pydantic import BaseModel, Field
from datetime import datetime

class CmdMetrics(BaseModel):
    cmd: Literal["extract","normalize","enrich","export","report","validate"]
    start_time: datetime = Field(default_factory=datetime.utcnow)
    end_time: Optional[datetime] = None
    duration_s: Optional[float] = None
    input: Dict[str, Any] = Field(default_factory=dict)
    counts: Dict[str, int] = Field(default_factory=dict)
    api: Dict[str, Any] = Field(default_factory=dict)
    cache: Dict[str, int] = Field(default_factory=dict)
    exit_code: int = 0
    notes: Optional[str] = None

    def finish(self, exit_code: int = 0, notes: Optional[str] = None):
        self.end_time = datetime.utcnow()
        if self.start_time and self.end_time:
            self.duration_s = (self.end_time - self.start_time).total_seconds()
        self.exit_code = exit_code
        if notes:
            self.notes = notes
        return self
