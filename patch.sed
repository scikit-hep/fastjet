s/from _fastjet import FastJetError/if __package__ or "." in __name__:\n    from ._fastjet import FastJetError\n  else:\n    from _fastjet import FastJetError/g
