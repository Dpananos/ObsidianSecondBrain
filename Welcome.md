Sounds like the experiment is tacking user activation and retention.  These sound like good primary metrics, but we need a better definition for these.  

Assuming I get to define what activation and retention means, I would say something like:

- Activation: % of users who make an edit to the template
- Retention: % of users who create a second file within some amount of time.
  
  
The key here is to create a metric which is sensitive enough but is not too noisy, irrelevant, or too far along in the user experience so as to be difficult to move.

Now, I need a single primary meric, so I am going to choose actication as my primary metric 


---

Maybe we need to refine activation so be edit and/or add to the template.  I suppose this depends on exactly how the events are set up, but generally we want to measure some sort of "usage" of the template without downright abandoning it.

Great point about retention -- maybe we need to refine this to be coming back to the dociument and making firther changes/edits?  Bascially something to measure if the user is getting value out of the document and comes back to reference/usie it.

I'm picking activation over retention because we can potentially get more signal here.  I'm assuming that retention is going to be smaller than activation and if both are noisy metrics then I would imagine retenion is going to require a larger sample size, meaning the experiment should run longer.  

Guard rails is a good point.  We might want to make sure users are not abandoning these docuiments, and we also want to make sure the product is smooth (so measuring something like latency is a good idea here.).  Thigns like crash rate and/or customer service metrics (e.g. reaching out to support) are not moved in he wrong direction

---

I'd probably randomize on the user here.  We're trying to move user level metrics (i.e. retention and activation).  .

Experiment duration can be handled with an MDE calculation.  We want to choose an MDE that is of a practical size, but is not so small that we are running the experiment too long.  If we have past experiments which tried to move activation at a similar point in the user flow, maybe we could do some clever math here but I think to a first approximation it is a good habit to specify a desired run time (e.g. 2 weeks), see what kind of MDEs we can get at that run time, then have a discussion on if we need more time and if that trade off is worth it or if we need to re-think our approach to decision making (e.g. maybe an experiment isn't feasible b/c of time constraints, what else can we do to get som feedback on if this is creating value).

So far as design pitfalls, I think one major thing I can think of is excluding users who already activated.  We maybe want to target new signups and randomize them at the right moment to measure if we are creating value for new users as opposed to diluting our experiment with users who couldn't get the outcome

---

Running experiments in increments of weeks allows you to marginalize (or, average over) day of week effects.  If you run only 5 od the 7 days, maybe you don't run on the weekend when activation is naturally lower and so the experiment looks overly good.

Good point on the control.  We may want to show some templates to users, but have them be randomly generated or sampled.  The idea is to try and determine if if is the curation of the templates as opposed to having templates (i.e. avoiding the blank slate syndrome)

---

Interesting results.  Activation is statistically significant meaning the treatment increased the conversion rate.

However, we are seeing that some of the guardrails are stat sig as well.  Looks like the curation lead to fewer templates browsed (maybe that is a good thing and we are surfacing the right templates).

Median page load time has increased.  This could signal that there is a bug or at least an inefficient implementation.  Maybe we should investigate.  In teh event we don't find osmething, we need to weigh the trade offs between higher activation rate and a slower product.

I'm unsurprised 7 day retention is not stat sig.  We said this was a low sensitivity metric and so this could be because we had insufficient power to detect an effect (though this does not mean there was NO effect, just that if ther ewas then we couldn't detect it)

---
The lift in activation is modest, yes.  Generally, filtering on stat sig results is subject to the winner's curse so while we have a stat sig result, maybe we should temper our expectations.

To answer your question, would I be excited on shipping this feaure if the lift was 0.8%?  probably not, that is a modest lift, and if that were the true lift then trh tradeoffs (esp with respect to latency) feel more clear cut

---

EXperimentation is a team sport.  I think if I had a PM who pushed back, it becomes a conversation more than a stats lesson.  There is real risk to shipping without further investigation into the latency bug. My recommendation is to fix the bug and re-test to validate we are seeing lift.

But at the end of the day, I am not here to enforce statsitics rules, just clarify what the risks are and what safe routes we can take.  

