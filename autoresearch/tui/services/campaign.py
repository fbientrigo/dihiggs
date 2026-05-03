"""
Campaign process service - manages supervisor subprocess lifecycle.
"""
import subprocess
import signal
import sys
from pathlib import Path
from typing import Optional, Callable
from dataclasses import dataclass, field
from enum import Enum


class ProcessState(Enum):
    IDLE = "idle"
    STARTING = "starting"
    RUNNING = "running"
    STOPPING = "stopping"
    STOPPED = "stopped"
    FAILED = "failed"


@dataclass
class CampaignProcess:
    """Manages a campaign supervisor subprocess."""
    
    config_path: Path
    outdir: Path
    state: ProcessState = ProcessState.IDLE
    process: Optional[subprocess.Popen] = None
    pid: Optional[int] = None
    error: Optional[str] = None
    _on_state_change: Optional[Callable[["CampaignProcess"], None]] = field(default=None, repr=False)
    
    def start(self) -> bool:
        """Start the campaign supervisor subprocess."""
        if self.state == ProcessState.RUNNING:
            return False
        
        if not self.config_path.exists():
            self.error = f"Config not found: {self.config_path}"
            self._set_state(ProcessState.FAILED)
            return False
        
        self._set_state(ProcessState.STARTING)
        
        try:
            # Launch supervisor as subprocess
            cmd = [
                sys.executable,
                "-m", "autoresearch.run_supervisor",
                str(self.config_path),
                "--outdir", str(self.outdir),
            ]
            
            self.process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                cwd=Path(__file__).parent.parent.parent.parent,  # repo root
            )
            self.pid = self.process.pid
            self._set_state(ProcessState.RUNNING)
            return True
            
        except Exception as e:
            self.error = str(e)
            self._set_state(ProcessState.FAILED)
            return False
    
    def stop(self, timeout: float = 10.0) -> bool:
        """Stop the campaign subprocess gracefully."""
        if self.process is None:
            self._set_state(ProcessState.STOPPED)
            return True
        
        self._set_state(ProcessState.STOPPING)
        
        try:
            # Send SIGTERM for graceful shutdown
            self.process.terminate()
            try:
                self.process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                # Force kill if graceful shutdown fails
                self.process.kill()
                self.process.wait(timeout=5.0)
            
            self._set_state(ProcessState.STOPPED)
            return True
            
        except Exception as e:
            self.error = str(e)
            self._set_state(ProcessState.FAILED)
            return False
    
    def poll(self) -> Optional[int]:
        """Check if process is still running. Returns exit code if finished."""
        if self.process is None:
            return None
        
        retcode = self.process.poll()
        if retcode is not None:
            if retcode == 0:
                self._set_state(ProcessState.STOPPED)
            else:
                self.error = f"Process exited with code {retcode}"
                self._set_state(ProcessState.FAILED)
        
        return retcode
    
    def is_running(self) -> bool:
        """Check if process is currently running."""
        if self.process is None:
            return False
        return self.process.poll() is None
    
    def _set_state(self, new_state: ProcessState):
        """Update state and notify callback."""
        self.state = new_state
        if self._on_state_change:
            self._on_state_change(self)


class CampaignService:
    """Service layer for managing campaign processes."""
    
    def __init__(self):
        self.active_campaign: Optional[CampaignProcess] = None
    
    def create_campaign(
        self,
        config_path: Path,
        outdir: Path,
        on_state_change: Optional[Callable[[CampaignProcess], None]] = None,
    ) -> CampaignProcess:
        """Create a new campaign process (does not start it)."""
        campaign = CampaignProcess(
            config_path=config_path,
            outdir=outdir,
            _on_state_change=on_state_change,
        )
        self.active_campaign = campaign
        return campaign
    
    def start_campaign(self) -> bool:
        """Start the active campaign."""
        if self.active_campaign is None:
            return False
        return self.active_campaign.start()
    
    def stop_campaign(self) -> bool:
        """Stop the active campaign."""
        if self.active_campaign is None:
            return False
        return self.active_campaign.stop()
    
    def get_state(self) -> ProcessState:
        """Get the current campaign state."""
        if self.active_campaign is None:
            return ProcessState.IDLE
        return self.active_campaign.state
